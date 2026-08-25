#!/usr/bin/env python3

import argparse
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract quality metrics from raw trio de novo candidates."
    )

    parser.add_argument(
        "--vcf",
        required=True,
        help="Input VCF containing raw de novo candidates"
    )

    parser.add_argument(
        "--ped", 
        required=True, 
        help="PED file defining the trio")

    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV containing candidate quality metrics"
    )

    return parser.parse_args()


def parse_info(info_string):
    """
    Convert the INFO field into a dictionary.
    Example:
        QD=10.2;MQ=59.4;FS=2.1
    becomes:
        {"QD": "10.2", "MQ": "59.4", "FS": "2.1"}
    """

    info = {}

    for item in info_string.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            # Flag-style INFO field, e.g. DB
            info[item] = True

    return info


def parse_sample(format_keys, sample_string):
    """
    Parse one sample FORMAT field into a dictionary.

    Example:
        FORMAT = GT:AD:DP:GQ:PL
        sample = 0/1:18,16:34:99:...
    """

    values = sample_string.split(":")

    return dict(zip(format_keys, values))


def safe_int(value):
    """
    Convert a value to integer when possible.
    Missing values become None.
    """

    if value in {None, ".", ""}:
        return None

    try:
        return int(value)
    except ValueError:
        return None


def safe_float(value):
    """
    Convert a value to float when possible.
    Missing values become None.
    """

    if value in {None, ".", ""}:
        return None

    try:
        return float(value)
    except ValueError:
        return None


def parse_ad(ad_value):
    """
    Parse biallelic allele depth.

    AD format:
        REF,ALT

    Returns:
        ref_depth, alt_depth
    """

    if ad_value in {None, ".", ""}:
        return None, None

    parts = ad_value.split(",")

    if len(parts) < 2:
        return None, None

    ref_depth = safe_int(parts[0])
    alt_depth = safe_int(parts[1])

    return ref_depth, alt_depth


def calculate_allele_balance(ref_depth, alt_depth):
    """
    Calculate ALT allele balance:

        ALT / (REF + ALT)
    """

    if ref_depth is None or alt_depth is None:
        return None

    total = ref_depth + alt_depth

    if total == 0:
        return None

    return alt_depth / total


def profile_candidates(vcf_path, ped_path,output_path):

    mother_index = None
    father_index = None
    proband_index = None

    rows_written = 0

    output_columns = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",

        "M_GT",
        "M_DP",
        "M_GQ",
        "M_REF",
        "M_ALT",

        "F_GT",
        "F_DP",
        "F_GQ",
        "F_REF",
        "F_ALT",

        "P_GT",
        "P_DP",
        "P_GQ",
        "P_REF",
        "P_ALT",
        "P_AB",

        "QD",
        "MQ",
        "FS",
        "SOR",
        "MQRankSum",
        "ReadPosRankSum"
    ]

    with open(vcf_path, "r") as vcf_file, open(
        output_path,
        "w",
        newline=""
    ) as output_file:

        writer = csv.DictWriter(
            output_file,
            fieldnames=output_columns,
            delimiter="\t"
        )

        writer.writeheader()

        for line in vcf_file:

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):

                header = line.rstrip("\n").split("\t")
                samples = header[9:]

                # These are the sample IDs from this specific trio.
                mother_id = father_id = proband_id = None
                for ped_line in open(ped_path):
                    f = ped_line.split()
                    if len(f) >= 6 and f[5] == "2":
                        proband_id, father_id, mother_id = f[1], f[2], f[3]
                if proband_id is None:
                    raise ValueError("No affected proband (phenotype=2) in PED.")
                

                if mother_id not in samples:
                    raise ValueError(f"{mother_id} not found in VCF.")

                if father_id not in samples:
                    raise ValueError(f"{father_id} not found in VCF.")

                if proband_id not in samples:
                    raise ValueError(f"{proband_id} not found in VCF.")

                mother_index = samples.index(mother_id) + 9
                father_index = samples.index(father_id) + 9
                proband_index = samples.index(proband_id) + 9

                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filter_status = fields[6]

            info = parse_info(fields[7])

            format_keys = fields[8].split(":")

            mother = parse_sample(
                format_keys,
                fields[mother_index]
            )

            father = parse_sample(
                format_keys,
                fields[father_index]
            )

            proband = parse_sample(
                format_keys,
                fields[proband_index]
            )

            # Parse allele depths
            m_ref, m_alt = parse_ad(mother.get("AD"))
            f_ref, f_alt = parse_ad(father.get("AD"))
            p_ref, p_alt = parse_ad(proband.get("AD"))

            # Proband heterozygous allele balance
            p_ab = calculate_allele_balance(
                p_ref,
                p_alt
            )

            row = {
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "QUAL": safe_float(qual),
                "FILTER": filter_status,

                "M_GT": mother.get("GT"),
                "M_DP": safe_int(mother.get("DP")),
                "M_GQ": safe_int(mother.get("GQ")),
                "M_REF": m_ref,
                "M_ALT": m_alt,

                "F_GT": father.get("GT"),
                "F_DP": safe_int(father.get("DP")),
                "F_GQ": safe_int(father.get("GQ")),
                "F_REF": f_ref,
                "F_ALT": f_alt,

                "P_GT": proband.get("GT"),
                "P_DP": safe_int(proband.get("DP")),
                "P_GQ": safe_int(proband.get("GQ")),
                "P_REF": p_ref,
                "P_ALT": p_alt,
                "P_AB": (
                    round(p_ab, 4)
                    if p_ab is not None
                    else None
                ),

                "QD": safe_float(info.get("QD")),
                "MQ": safe_float(info.get("MQ")),
                "FS": safe_float(info.get("FS")),
                "SOR": safe_float(info.get("SOR")),
                "MQRankSum": safe_float(info.get("MQRankSum")),
                "ReadPosRankSum": safe_float(
                    info.get("ReadPosRankSum")
                )
            }

            writer.writerow(row)
            rows_written += 1

    print(f"Candidates profiled: {rows_written}")
    print(f"Output written to: {output_path}")


def main():

    args = parse_args()

    profile_candidates(
        args.vcf,
        args.ped,
        args.output
    )


if __name__ == "__main__":
    main()