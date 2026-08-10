#!/usr/bin/env bash

set -euo pipefail

VCF=${VCF:-data/biallelic_CDL-068-99.vcf}
PED=${PED:-data/CDL-068-99.ped}
OUT_DIR=${OUT_DIR:-results}

mkdir -p "$OUT_DIR"

echo "=== Step 1: Identify strict de novo candidates ==="
python3 scripts/find_denovo_candidates.py \
  --vcf "$VCF" \
  --ped "$PED" \
  --output "$OUT_DIR/denovo_candidates.vcf"

echo "=== Step 2: Genotype-level QC ==="
python3 scripts/filter_genotype_qc.py \
  --vcf "$VCF" \
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