#!/usr/bin/env python3

import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply genotype-level QC to trio de novo SNV candidates."
    )

    parser.add_argument(
        "--vcf",
        required=True,
        help="Input biallelic VCF"
    )

    parser.add_argument(
        "--ped",
        required=True,
        help="Input PED file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output VCF containing genotype-QC-passing de novo SNVs"
    )

    parser.add_argument(
        "--summary",
        required=True,
        help="Output text summary of filtering counts"
    )

    return parser.parse_args()


def read_pedigree(ped_path):
    """
    Identify the affected proband (phenotype=2) and parents.
    """

    affected = []

    with open(ped_path) as ped_file:
        for line in ped_file:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 6:
                raise ValueError("PED must contain at least 6 columns.")

            family, sample, father, mother, sex, phenotype = fields[:6]

            if phenotype == "2":
                affected.append(
                    {
                        "proband": sample,
                        "father": father,
                        "mother": mother
                    }
                )

    if len(affected) != 1:
        raise ValueError(
            "Expected exactly one affected proband (phenotype=2)."
        )

    return affected[0]


def normalize_gt(gt):
    if gt in {None, ".", "./.", ".|."}:
        return None

    gt = gt.replace("|", "/")

    alleles = gt.split("/")

    if len(alleles) != 2:
        return gt

    return "/".join(sorted(alleles))


def parse_sample(format_keys, sample_string):
    values = sample_string.split(":")
    return dict(zip(format_keys, values))


def parse_int(value):
    if value in {None, ".", ""}:
        return None

    try:
        return int(value)
    except ValueError:
        return None


def parse_ad(ad):
    """
    Biallelic AD:
        REF,ALT
    """

    if ad in {None, ".", ""}:
        return None, None

    values = ad.split(",")

    if len(values) < 2:
        return None, None

    return parse_int(values[0]), parse_int(values[1])


def calc_vaf(ref_depth, alt_depth):
    if ref_depth is None or alt_depth is None:
        return None

    total = ref_depth + alt_depth

    if total == 0:
        return None

    return alt_depth / total


def pl_supports(pl_string, expected_index):
    """
    Biallelic diploid PL order:

        PL[0] = 0/0
        PL[1] = 0/1
        PL[2] = 1/1
    """

    if pl_string in {None, ".", ""}:
        return False

    try:
        values = [int(x) for x in pl_string.split(",")]
    except ValueError:
        return False

    if len(values) < 3:
        return False

    return values[expected_index] == min(values[:3])


def main():
    args = parse_args()

    trio = read_pedigree(args.ped)

    mother_id = trio["mother"]
    father_id = trio["father"]
    proband_id = trio["proband"]

    print(f"Mother:  {mother_id}")
    print(f"Father:  {father_id}")
    print(f"Proband: {proband_id}")

    sample_indices = {}

    counts = {
        "strict_denovo": 0,
        "ordinary_snv": 0,
        "depth": 0,
        "gq": 0,
        "proband_vaf": 0,
        "parent_vaf": 0,
        "pl": 0
    }

    with open(args.vcf) as vcf_file, open(args.output, "w") as out_vcf:

        for line in vcf_file:

            # Preserve metadata
            if line.startswith("##"):
                out_vcf.write(line)
                continue

            # Parse sample order
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]

                for role, sample_id in {
                    "mother": mother_id,
                    "father": father_id,
                    "proband": proband_id
                }.items():

                    if sample_id not in samples:
                        raise ValueError(
                            f"{role} sample '{sample_id}' not found in VCF."
                        )

                    sample_indices[role] = samples.index(sample_id) + 9

                out_vcf.write(line)
                continue

            if line.startswith("#"):
                out_vcf.write(line)
                continue

            fields = line.rstrip("\n").split("\t")

            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]

            format_keys = fields[8].split(":")

            mother = parse_sample(
                format_keys,
                fields[sample_indices["mother"]]
            )

            father = parse_sample(
                format_keys,
                fields[sample_indices["father"]]
            )

            proband = parse_sample(
                format_keys,
                fields[sample_indices["proband"]]
            )

            m_gt = normalize_gt(mother.get("GT"))
            f_gt = normalize_gt(father.get("GT"))
            p_gt = normalize_gt(proband.get("GT"))

            # -------------------------------------------------
            # 1. Strict de novo genotype
            # -------------------------------------------------

            if not (
                m_gt == "0/0"
                and f_gt == "0/0"
                and p_gt == "0/1"
            ):
                continue

            counts["strict_denovo"] += 1

            # -------------------------------------------------
            # 2. Ordinary SNVs only
            # -------------------------------------------------

            if alt == "*":
                continue

            if len(ref) != 1 or len(alt) != 1:
                continue

            counts["ordinary_snv"] += 1

            # Extract DP/GQ
            m_dp = parse_int(mother.get("DP"))
            f_dp = parse_int(father.get("DP"))
            p_dp = parse_int(proband.get("DP"))

            m_gq = parse_int(mother.get("GQ"))
            f_gq = parse_int(father.get("GQ"))
            p_gq = parse_int(proband.get("GQ"))

            # Extract AD
            m_ref, m_alt = parse_ad(mother.get("AD"))
            f_ref, f_alt = parse_ad(father.get("AD"))
            p_ref, p_alt = parse_ad(proband.get("AD"))

            m_vaf = calc_vaf(m_ref, m_alt)
            f_vaf = calc_vaf(f_ref, f_alt)
            p_vaf = calc_vaf(p_ref, p_alt)

            # -------------------------------------------------
            # 3. DP >= 10 in trio
            # -------------------------------------------------

            if None in {m_dp, f_dp, p_dp}:
                continue

            if not (
                m_dp >= 10
                and f_dp >= 10
                and p_dp >= 10
            ):
                continue

            counts["depth"] += 1

            # -------------------------------------------------
            # 4. GQ >= 30 in trio
            # -------------------------------------------------

            if None in {m_gq, f_gq, p_gq}:
                continue

            if not (
                m_gq >= 30
                and f_gq >= 30
                and p_gq >= 30
            ):
                continue

            counts["gq"] += 1

            # -------------------------------------------------
            # 5. Proband VAF 0.30-0.70
            # -------------------------------------------------

            if p_vaf is None:
                continue

            if not (0.30 <= p_vaf <= 0.70):
                continue

            counts["proband_vaf"] += 1

            # -------------------------------------------------
            # 6. Parent ALT VAF < 0.10
            # -------------------------------------------------

            if m_vaf is None or f_vaf is None:
                continue

            if not (
                m_vaf < 0.10
                and f_vaf < 0.10
            ):
                continue

            counts["parent_vaf"] += 1

            # -------------------------------------------------
            # 7. PL consistency
            # -------------------------------------------------

            if not (
                pl_supports(mother.get("PL"), 0)
                and pl_supports(father.get("PL"), 0)
                and pl_supports(proband.get("PL"), 1)
            ):
                continue

            counts["pl"] += 1

            # Candidate passes genotype-level QC
            out_vcf.write(line)

    summary_lines = [
        ("Initial strict genotype candidates", counts["strict_denovo"]),
        ("Ordinary SNV candidates", counts["ordinary_snv"]),
        ("DP >= 10 in trio", counts["depth"]),
        ("GQ >= 30 in trio", counts["gq"]),
        ("Proband VAF 0.30-0.70", counts["proband_vaf"]),
        ("Parental ALT fraction < 0.10", counts["parent_vaf"]),
        ("PL supports expected genotype", counts["pl"])
    ]

    print("\n=== GENOTYPE-LEVEL QC ===\n")

    for label, count in summary_lines:
        print(f"{label:<40} {count}")

    with open(args.summary, "w") as summary_file:
        summary_file.write(
            "Genotype-level de novo QC summary\n"
        )
        summary_file.write(
            "=================================\n\n"
        )

        for label, count in summary_lines:
            summary_file.write(
                f"{label:<40} {count}\n"
            )


if __name__ == "__main__":
    main()