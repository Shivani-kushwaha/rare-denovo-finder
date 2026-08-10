# De novo variant identification

- A reproducible workflow for identifying **strict de novo SNV candidates** from a trio whole-genome VCF and filtering them to a high-confidence candidate set using genotype-level, site-level, and cohort-level evidence.
- The analysis uses the trio probing, father, and mother. The supplied VCF was decomposed into biallelic records using bcftools norm -m -.

## Method summary
**Task 1 — Strict de novo identification.**  
Variants are selected when both parents are homozygous reference (0/0) and the proband is heterozygous (0/1). Symbolic spanning-deletion alleles (*) are excluded and the analysis is restricted to SNVs.

**Task 2 — Genotype-level QC.**  
Candidates require Genotype-level filters
- GT: parents 0/0, proband 0/1
- DP >= 10 in all trio members
- GQ >= 30
- Proband VAF 0.30–0.70
- Parent ALT VAF < 0.10
- PL supports expected genotype

**Task 2 — Site-level QC.**  
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

**Task 2 — De novo prioritization.**  
After technical QC, cohort singletons (AC=1) are prioritized as the most likely private de novo events. This is used as a prioritization criterion rather than an absolute definition of a de novo mutation.

## Quick start
- Place the input files in the data folder:
1. data/biallelic.vcf
2. data/sample.ped

- Run the complete workflow from the repository root:
bash scripts/run_all.sh

The wrapper executes:
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

Individual scripts can also be run separately if a specific filtering stage needs to be inspected.

## Results
Filtering step                                 Remaining Variant 
Strict Mendelian-pattern candidates               446
SNVs after excluding symbolic '*' alleles         210
Genotype-QC passing                               179
Site-QC passing                                   169
Prioritized cohort singletons (AC=1)               90

## Repository structure
rare_denovo-finder
├── README.md
├── data/
│   ├── biallelic.vcf
│   └── sample.ped
├── scripts/
│   ├── find_denovo_candidates.py
│   ├── filter_genotype_qc.py
│   ├── filter_site_qc.py
│   ├── filter_denovo_singletons.py
│   └── run_all.sh
└── results/
    ├── denovo_candidates.vcf
    ├── genotype_qc_candidates.vcf
    ├── site_qc_candidates.vcf
    └── final_denovo_candidates.vcf

- run_all.sh reproduces the final filtering workflow; metric profiling scripts are provided separately to inspect distributions used to justify thresholds.
- Tasks 3–6, including the expected number of true-positive de novo variants, sources of false-positive calls, expected structural variant burden, and CNV filtering strategy, are discussed in the accompanying report.
