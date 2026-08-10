#!/usr/bin/env python3

import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply site-level QC filters to genotype-QC de novo SNVs."
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Input VCF containing genotype-QC-passing de novo SNVs"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output VCF containing variants passing site-level QC"
    )
    parser.add_argument(
        "--summary",
        required=True,
        help="Output text file containing site-level filtering counts"
    )
    return parser.parse_args()


def parse_info(info_string):
    """
    Parse VCF INFO field into a dictionary.
    """

    info = {}
    for item in info_string.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True

    return info


def safe_float(value):
    """
    Convert a VCF value to float.
    Return None for missing/unparseable values.
    """

    if value in {None, ".", ""}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def main():
    args = parse_args()

    counts = {
        "input": 0,
        "pass": 0,
        "qual": 0,
        "qd": 0,
        "mq": 0,
        "fs": 0,
        "sor": 0,
        "mq_rank_sum": 0,
        "read_pos_rank_sum": 0
    }

    header_lines = []
    passing_lines = []

    with open(args.vcf) as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                header_lines.append(line)
                continue
            counts["input"] += 1

            fields = line.rstrip("\n").split("\t")

            filter_status = fields[6]
            qual = safe_float(fields[5])
            info = parse_info(fields[7])

            qd = safe_float(info.get("QD"))
            mq = safe_float(info.get("MQ"))
            fs = safe_float(info.get("FS"))
            sor = safe_float(info.get("SOR"))
            mq_rank_sum = safe_float(info.get("MQRankSum"))
            read_pos_rank_sum = safe_float(
                info.get("ReadPosRankSum")
            )

            # -------------------------------------------------
            # 1. Existing VCF FILTER status
            # -------------------------------------------------

            if filter_status != "PASS":
                continue
            counts["pass"] += 1

            # -------------------------------------------------
            # 2. QUAL >= 30
            # -------------------------------------------------

            if qual is None or qual < 30:
                continue
            counts["qual"] += 1

            # -------------------------------------------------
            # 3. QD >= 2
            # -------------------------------------------------

            if qd is None or qd < 2.0:
                continue
            counts["qd"] += 1

            # -------------------------------------------------
            # 4. MQ >= 40
            # -------------------------------------------------

            if mq is None or mq < 40.0:
                continue
            counts["mq"] += 1

            # -------------------------------------------------
            # 5. FS <= 60
            # -------------------------------------------------

            if fs is None or fs > 60.0:
                continue
            counts["fs"] += 1

            # -------------------------------------------------
            # 6. SOR <= 3
            # -------------------------------------------------

            if sor is None or sor > 3.0:
                continue
            counts["sor"] += 1

            # -------------------------------------------------
            # 7. MQRankSum >= -12.5
            # -------------------------------------------------

            if mq_rank_sum is None or mq_rank_sum < -12.5:
                continue
            counts["mq_rank_sum"] += 1

            # -------------------------------------------------
            # 8. ReadPosRankSum >= -8
            # -------------------------------------------------

            if (
                read_pos_rank_sum is None
                or read_pos_rank_sum < -8.0
            ):
                continue
            counts["read_pos_rank_sum"] += 1

            # Variant passes all site-level QC
            passing_lines.append(line)

    # ---------------------------------------------------------
    # Write final VCF
    # ---------------------------------------------------------

    with open(args.output, "w") as out_vcf:
        for line in header_lines:
            out_vcf.write(line)
        for line in passing_lines:
            out_vcf.write(line)

    # ---------------------------------------------------------
    # Create summary
    # ---------------------------------------------------------

    summary_lines = [
        ("Input genotype-QC candidates", counts["input"]),
        ("FILTER = PASS", counts["pass"]),
        ("QUAL >= 30", counts["qual"]),
        ("QD >= 2", counts["qd"]),
        ("MQ >= 40", counts["mq"]),
        ("FS <= 60", counts["fs"]),
        ("SOR <= 3", counts["sor"]),
        ("MQRankSum >= -12.5", counts["mq_rank_sum"]),
        ("ReadPosRankSum >= -8", counts["read_pos_rank_sum"])
    ]

    print("\n=== SITE-LEVEL QC ===\n")

    for label, count in summary_lines:
        print(f"{label:<36} {count}")
    final_count = len(passing_lines)
    print()
    print(f"Final site-QC candidates: {final_count}")

    # ---------------------------------------------------------
    # Write summary file
    # ---------------------------------------------------------

    with open(args.summary, "w") as summary_file:

        summary_file.write(
            "Site-level de novo QC summary\n"
        )
        summary_file.write(
            "============================\n\n"
        )
        for label, count in summary_lines:
            summary_file.write(
                f"{label:<36} {count}\n"
            )
        summary_file.write(
            f"\nFinal site-QC candidates: {final_count}\n"
        )


if __name__ == "__main__":
    main()