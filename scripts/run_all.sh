#!/usr/bin/env bash

set -euo pipefail

VCF=${VCF:-}
PED=${PED:-}
OUT_DIR=${OUT_DIR:-results}

if [[ -z "$VCF" || -z "$PED" ]]; then
  echo "Usage: VCF=<trio.vcf> PED=<trio.ped> [OUT_DIR=results] bash scripts/run_all.sh" >&2
  exit 1
fi

for f in "$VCF" "$PED"; do
  if [[ ! -f "$f" ]]; then
    echo "Error: input file not found: $f" >&2
    exit 1
  fi
done

mkdir -p "$OUT_DIR"

echo "=== Step 1: Identify strict de novo candidates ==="
python3 scripts/find_denovo_candidates.py \
  --vcf "$VCF" \
  --ped "$PED" \
  --output "$OUT_DIR/denovo_candidates.vcf"

echo "=== Step 2: Genotype-level QC ==="
python3 scripts/filter_genotype_qc.py \
  --vcf "$OUT_DIR/denovo_candidates.vcf" \
  --ped "$PED" \
  --output "$OUT_DIR/genotype_qc_candidates.vcf" \
  --summary "$OUT_DIR/genotype_qc_summary.txt"

echo "=== Step 3: Site-level QC ==="
python3 scripts/filter_site_qc.py \
  --vcf "$OUT_DIR/genotype_qc_candidates.vcf" \
  --output "$OUT_DIR/site_qc_candidates.vcf" \
  --summary "$OUT_DIR/site_qc_summary.txt"

echo "=== Step 4: De novo singleton prioritization ==="
python3 scripts/filter_denovo_singletons.py \
  --vcf "$OUT_DIR/site_qc_candidates.vcf" \
  --output "$OUT_DIR/final_denovo_candidates.vcf" \
  --summary "$OUT_DIR/final_denovo_summary.txt"

echo
echo "Analysis complete."
echo "Final candidates: $OUT_DIR/final_denovo_candidates.vcf"