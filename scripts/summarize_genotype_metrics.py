#!/usr/bin/env python3

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Summarize genotype-level metrics for de novo SNV candidates."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input denovo_metrics.tsv file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV summary file"
    )

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    # Remove spanning-deletion symbolic alleles
    df = df[df["ALT"] != "*"].copy()

    # Calculate parental ALT allele fractions
    df["M_VAF"] = df["M_ALT"] / (df["M_REF"] + df["M_ALT"])
    df["F_VAF"] = df["F_ALT"] / (df["F_REF"] + df["F_ALT"])

    stats_order = ["min", "25%", "50%", "75%", "max"]

    rows = []

    # Candidate count
    rows.append({
        "SECTION": "CANDIDATE_COUNT",
        "METRIC": "SNV_candidates_analyzed",
        "STAT": "count",
        "MOTHER": "",
        "FATHER": "",
        "PROBAND": len(df)
    })

    # Depth
    dp = (
        df[["M_DP", "F_DP", "P_DP"]]
        .describe()
        .loc[stats_order]
    )

    for stat in stats_order:
        rows.append({
            "SECTION": "DEPTH_DP",
            "METRIC": "DP",
            "STAT": stat,
            "MOTHER": dp.loc[stat, "M_DP"],
            "FATHER": dp.loc[stat, "F_DP"],
            "PROBAND": dp.loc[stat, "P_DP"]
        })

    # Genotype quality
    gq = (
        df[["M_GQ", "F_GQ", "P_GQ"]]
        .describe()
        .loc[stats_order]
    )

    for stat in stats_order:
        rows.append({
            "SECTION": "GENOTYPE_QUALITY_GQ",
            "METRIC": "GQ",
            "STAT": stat,
            "MOTHER": gq.loc[stat, "M_GQ"],
            "FATHER": gq.loc[stat, "F_GQ"],
            "PROBAND": gq.loc[stat, "P_GQ"]
        })

    # Proband ALT allele fraction
    p_vaf = (
        df["P_AB"]
        .describe()
        .loc[stats_order]
    )

    for stat in stats_order:
        rows.append({
            "SECTION": "PROBAND_ALT_ALLELE_FRACTION",
            "METRIC": "P_VAF",
            "STAT": stat,
            "MOTHER": "",
            "FATHER": "",
            "PROBAND": p_vaf.loc[stat]
        })

    # Parent ALT allele fraction
    parent_vaf = (
        df[["M_VAF", "F_VAF"]]
        .describe()
        .loc[stats_order]
    )

    for stat in stats_order:
        rows.append({
            "SECTION": "PARENT_ALT_ALLELE_FRACTION",
            "METRIC": "VAF",
            "STAT": stat,
            "MOTHER": parent_vaf.loc[stat, "M_VAF"],
            "FATHER": parent_vaf.loc[stat, "F_VAF"],
            "PROBAND": ""
        })

    # Count candidates with any parental ALT reads
    parental_alt_count = (
        (df["M_ALT"] > 0) |
        (df["F_ALT"] > 0)
    ).sum()

    rows.append({
        "SECTION": "PARENT_ALT_SUPPORT",
        "METRIC": "Candidates_with_any_parental_ALT_reads",
        "STAT": "count",
        "MOTHER": "",
        "FATHER": "",
        "PROBAND": parental_alt_count
    })

    summary_df = pd.DataFrame(rows)

    summary_df.to_csv(
        args.output,
        sep="\t",
        index=False
    )

    print(f"SNV candidates analyzed: {len(df)}")
    print(f"Summary written to: {args.output}")


if __name__ == "__main__":
    main()