Genotype-level filters
- GT: parents 0/0, proband 0/1
- DP >= 10 in all trio members
- GQ >= 30
- Proband VAF 0.30–0.70
- Parent ALT VAF < 0.10
- PL supports expected genotype

Site-level filters
- FILTER = PASS
- QUAL >= 30
- QD >= 2
- MQ >= 40
- FS <= 60
- SOR <= 3
- MQRankSum >= -12.5
- ReadPosRankSum >= -8


SUMMARY:
Strict de novo candidates: 446
Ordinary SNVs: 210
After genotype-level QC: 179
After site-level QC: 169