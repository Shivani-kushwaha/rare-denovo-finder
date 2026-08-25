# De Novo Variant Identification

- A reproducible workflow for identifying **strict de novo SNV candidates** from a trio whole-genome VCF and filtering them to a high-confidence candidate set using genotype-level, site-level, and cohort-level evidence.
- The analysis uses the trio proband, father, and mother. The supplied VCF was decomposed into biallelic records using bcftools norm -m -.

## Method summary
**Task 1 — Strict de novo identification.**  
Variants are selected when both parents are homozygous reference (0/0) and the proband is heterozygous (0/1). Symbolic spanning-deletion alleles (*) are excluded and the analysis is restricted to SNVs.

**Task 2a — Genotype-level QC.**  
Candidates require Genotype-level filters
- GT: parents 0/0, proband 0/1
- DP >= 10 in all trio members
- GQ >= 30
- Proband VAF 0.30–0.70
- Parent ALT VAF < 0.10
- PL supports expected genotype

**Task 2b — Site-level QC.**  
Variants require site-level filters
- FILTER = PASS
- QUAL >= 30
- QD >= 2
- MQ >= 40
- FS <= 60
- SOR <= 3
- MQRankSum >= -12.5
- ReadPosRankSum >= -8
Because the supplied VCF had already undergone VQSR, no additional arbitrary VQSLOD threshold was applied.

**Task 2c — De novo prioritization.**  
After technical QC, cohort singletons (AC=1) are prioritized as the most likely private de novo events. This is used as a prioritization criterion rather than an absolute definition of a de novo mutation.

## Quick start
Place the input files in the 'data/' folder:
- 'data/biallelic_CDL-068-99.vcf'
- 'data/CDL-068-99.ped'

Run the complete workflow from the repository root:
```bash
bash scripts/run_all.sh
```

To run on a different trio, override the defaults:
```bash
VCF=path/to/your.vcf PED=path/to/your.ped bash scripts/run_all.sh
```
  
## The wrapper executes:
```
Input VCF
   ↓
Strict de novo identification
   ↓
Genotype-level QC
   ↓
Site-level QC
   ↓
Singleton (AC=1) prioritization
   ↓
Final candidate VCF
```

## Results
Filtering step                                 Remaining Variant 
Strict Mendelian-pattern candidates               446
SNVs after excluding symbolic '*' alleles         210
Genotype-QC passing                               179
Site-QC passing                                   169
**Prioritized cohort singletons (AC=1)            90**

The 79 site-QC candidates with AC > 1 were removed at the final step (169 - 79 = 90).

## Repository structure

```
rare-denovo-finder/
├── README.md
├── environment.yml
├── LICENSE
├── .gitignore
├── data/
│   ├── biallelic_CDL-068-99.vcf
│   └── CDL-068-99.ped
├── docs/
│   └── filtering_rationale.md
├── scripts/
│   ├── find_denovo_candidates.py      Extracts the strict de novo genotype pattern (parents 0/0, proband 0/1)
│   ├── filter_genotype_qc.py          Applies trio genotype-level QC: depth, GQ, proband VAF, parental ALT fraction, PL
│   ├── filter_site_qc.py              Applies GATK best-practice site-level hard filters
│   ├── filter_denovo_singletons.py    Prioritizes cohort singletons (AC=1) to produce the final candidate set
│   ├── run_all.sh                     Runs the four filtering stages in order
│   ├── profile_denovo_metrics.py      Extracts per-variant quality metrics to TSV for threshold selection
│   ├── summarize_genotype_metrics.py  Reports genotype-metric distributions from the metrics TSV
│   └── summarize_site_metrics.py      Reports site-annotation distributions across the genotype-QC set
└── results/                           
```

## Notes
- `run_all.sh` reproduces the filtering workflow. The profiling and summary scripts are
  run separately to inspect the metric distributions used to justify each threshold,
  rather than adopting published cutoffs unexamined.
- Tasks 3–6 — expected number of true-positive de novo variants, sources of false-positive
  calls, expected structural variant burden, and CNV filtering strategy — are discussed in
  the accompanying report.
