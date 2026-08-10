#!/usr/bin/env python3

import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Retain singleton de novo SNV candidates (AC=1) "
            "after genotype- and site-level QC."
        )
    )

    parser.add_argument(
        "--vcf",
        required=True,
        help="Input site-QC-passing VCF"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output VCF containing singleton de novo candidates"
    )

    parser.add_argument(
        "--summary",
        required=True,
        help="Output text summary"
    )

    return parser.parse_args()


def parse_info(info_string):
    info = {}

    for item in info_string.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True

    return info


def main():
    args = parse_args()

    total = 0
    singleton_count = 0
    non_singleton_count = 0

    with open(args.vcf) as in_vcf, open(args.output, "w") as out_vcf:

        for line in in_vcf:

            if line.startswith("#"):
                out_vcf.write(line)
                continue

            total += 1

            fields = line.rstrip("\n").split("\t")
            info = parse_info(fields[7])

            ac = info.get("AC")

            if ac is None:
                raise ValueError(
                    "AC annotation is missing from one or more variants."
                )

            try:
                ac = int(ac)
            except ValueError:
                raise ValueError(
                    f"Unexpected AC value: {ac}"
                )

            if ac == 1:
                out_vcf.write(line)
                singleton_count += 1
            else:
                non_singleton_count += 1

    print("\n=== DE NOVO SINGLETON FILTER ===\n")
    print(f"Input site-QC candidates : {total}")
    print(f"AC = 1 singleton SNVs    : {singleton_count}")
    print(f"AC > 1 removed           : {non_singleton_count}")

    with open(args.summary, "w") as summary:

        summary.write("De novo singleton filtering summary\n")
        summary.write("===================================\n\n")

        summary.write(
            f"Input site-QC candidates : {total}\n"
        )
        summary.write(
            f"AC = 1 singleton SNVs    : {singleton_count}\n"
        )
        summary.write(
            f"AC > 1 removed           : {non_singleton_count}\n"
        )


if __name__ == "__main__":
    main()