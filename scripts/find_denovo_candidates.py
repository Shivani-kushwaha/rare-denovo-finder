#!/usr/bin/env python3
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Identify strict de novo variants from a biallelic trio VCF. "
            "A strict de novo candidate is defined as 0/0 in both parents and 0/1 in the proband "
        )
    )

    parser.add_argument(
        "--vcf",
        required=True,
        help="Input biallelic VCF file"
    )

    parser.add_argument(
        "--ped",
        required=True,
        help="Input PED file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output VCF containing strict de novo candidates"
    )

    return parser.parse_args()


def read_pedigree(ped_path):
    """
    Read PED file and identify the affected proband and their parents.
    The affected individual is treated as the proband.
    """

    affected_individuals = []

    with open(ped_path, "r") as ped_file:
        for line in ped_file:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 6:
                raise ValueError(
                    "PED file must contain at least 6 columns."
                )

            family_id, sample_id, father_id, mother_id, sex, phenotype = fields[:6]

            if phenotype == "2":
                affected_individuals.append(
                    {
                        "proband": sample_id,
                        "father": father_id,
                        "mother": mother_id
                    }
                )

    if len(affected_individuals) == 0:
        raise ValueError(
            "No affected individual (phenotype=2) found in PED file."
        )

    if len(affected_individuals) > 1:
        raise ValueError(
            "More than one affected individual found. "
            "Proband cannot be determined unambiguously."
        )

    trio = affected_individuals[0]

    if trio["father"] in {"0", "."} or trio["mother"] in {"0", "."}:
        raise ValueError(
            "Affected individual does not have both parents specified."
        )

    return trio


def normalize_genotype(gt):
    """
    Normalize phased genotypes to unphased representation.

    Examples:
        0|0 -> 0/0
        0|1 -> 0/1
        1|0 -> 0/1
    """

    if gt is None:
        return None

    gt = gt.replace("|", "/")

    if gt in {".", "./."}:
        return None

    alleles = gt.split("/")

    if len(alleles) != 2:
        return gt

    return "/".join(sorted(alleles))


def get_genotype(sample_field, format_keys):
    """
    Extract GT from a sample field using the FORMAT column.
    """

    values = sample_field.split(":")

    sample_data = dict(zip(format_keys, values))

    return normalize_genotype(sample_data.get("GT"))


def identify_denovo(vcf_path, ped_path, output_path):
    trio = read_pedigree(ped_path)

    mother_id = trio["mother"]
    father_id = trio["father"]
    proband_id = trio["proband"]

    print(f"Mother:  {mother_id}")
    print(f"Father:  {father_id}")
    print(f"Proband: {proband_id}")

    sample_indices = {}

    total_variants = 0
    denovo_count = 0

    with open(vcf_path, "r") as vcf_file, open(output_path, "w") as output_file:

        for line in vcf_file:

            # Preserve metadata lines
            if line.startswith("##"):
                output_file.write(line)
                continue

            # Parse VCF header
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")

                samples = header[9:]

                required_samples = {
                    "mother": mother_id,
                    "father": father_id,
                    "proband": proband_id
                }

                for role, sample_id in required_samples.items():
                    if sample_id not in samples:
                        raise ValueError(
                            f"{role.capitalize()} sample '{sample_id}' "
                            "was not found in the VCF."
                        )

                    # +9 because sample columns begin at VCF column 10
                    sample_indices[role] = samples.index(sample_id) + 9

                output_file.write(line)
                continue

            if line.startswith("#"):
                output_file.write(line)
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                continue

            total_variants += 1

            format_keys = fields[8].split(":")

            mother_gt = get_genotype(
                fields[sample_indices["mother"]],
                format_keys
            )

            father_gt = get_genotype(
                fields[sample_indices["father"]],
                format_keys
            )

            proband_gt = get_genotype(
                fields[sample_indices["proband"]],
                format_keys
            )

            # Strict de novo definition for a biallelic VCF
            if (
                mother_gt == "0/0"
                and father_gt == "0/0"
                and proband_gt == "0/1"
            ):
                output_file.write(line)
                denovo_count += 1

    print()
    print(f"Total variants examined: {total_variants}")
    print(f"Strict de novo candidates: {denovo_count}")
    print(f"Output written to: {output_path}")


def main():
    args = parse_args()

    identify_denovo(
        vcf_path=args.vcf,
        ped_path=args.ped,
        output_path=args.output
    )


if __name__ == "__main__":
    main()