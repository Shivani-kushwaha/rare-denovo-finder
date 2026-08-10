#!/usr/bin/env python3

import argparse


def parse_info(info_string):
    info = {}

    for item in info_string.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True

    return info


def safe_float(value):
    if value in {None, ".", ""}:
        return None

    try:
        return float(value)
    except ValueError:
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Summarize site-level metrics in genotype-QC de novo SNVs."
    )

    parser.add_argument("--vcf", required=True)

    args = parser.parse_args()

    metrics = {
        "QUAL": [],
        "QD": [],
        "MQ": [],
        "FS": [],
        "SOR": [],
        "BaseQRankSum": [],
        "MQRankSum": [],
        "ReadPosRankSum": [],
        "VQSLOD": []
    }

    n = 0

    with open(args.vcf) as f:

        for line in f:

            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            n += 1

            qual = safe_float(fields[5])

            if qual is not None:
                metrics["QUAL"].append(qual)

            info = parse_info(fields[7])

            for metric in metrics:

                if metric == "QUAL":
                    continue

                value = safe_float(info.get(metric))

                if value is not None:
                    metrics[metric].append(value)

    print(f"\nSite-QC candidates analyzed: {n}")

    print("\n=== SITE-LEVEL METRICS ===")

    for metric, values in metrics.items():

        print(f"\n{metric}")

        if not values:
            print("  No values available")
            continue

        values = sorted(values)

        def percentile(p):
            index = int((len(values) - 1) * p)
            return values[index]

        print(f"  min : {min(values):.3f}")
        print(f"  25% : {percentile(0.25):.3f}")
        print(f"  50% : {percentile(0.50):.3f}")
        print(f"  75% : {percentile(0.75):.3f}")
        print(f"  max : {max(values):.3f}")


if __name__ == "__main__":
    main()