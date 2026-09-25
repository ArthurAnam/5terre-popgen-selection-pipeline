# Methods overview

## Cinque Terre dataset

- WGS cohort
- 50 starting samples
- all males
- approximately 30x mean depth
- GRCh37.p13 / hs37d5
- delivered VCF already restricted upstream to PASS biallelic SNPs

Variant calling, joint genotyping, SNP selection, and VQSR were performed upstream.

Four reviewed sample-level outliers are excluded, leaving **46 samples**.

## CT variant QC

Core autosomal QC:

- chromosomes 1-22
- remove monomorphic sites
- remove A/T and C/G strand-ambiguous SNPs
- remove sites with missingness >5%
- remove individuals with missingness >5%
- repeat cohort-dependent site QC after any sample removal
- chromosome-wise HWE Bonferroni threshold: `0.05 / n_tested_chr`

No additional arbitrary QUAL, DP, GQ, or allele-balance hard filters are imposed after sample review.

Final CT core-QC panel:

- 46 samples
- 9,000,246 SNPs

## 1000 Genomes EUR

Phase 3 EUR reference:

- CEU: 99
- FIN: 99
- GBR: 91
- IBS: 107
- TSI: 107
- total: 503

Roles:

- CT + EUR PCA / ROH context
- EUR reference for SHAPEIT2
- CEU, TSI, IBS as production selection comparators
- FIN and GBR as LD-context populations only

For cross-dataset harmonization, only exact `CHR:POS:REF:ALT` matches are retained. A/T and C/G SNPs are removed before harmonization.

## Relatedness

KING v2.3.2 is run within CT.

No pair exceeds the 0.0442 third-degree kinship boundary; no sample is excluded on relatedness grounds.

## PCA

Joint CT + EUR PCA:

- samples: 549
- joint MAF >=0.05 panel: 4,878,327 SNPs
- LD-pruned panel: 295,375 SNPs
- high-LD-region sensitivity panel: 287,815 SNPs

The first three PCs are essentially unchanged by the high-LD-region sensitivity.
For PC1-PC2, a matched-sample Procrustes comparison of the masked and unmasked
solutions gives a standardized disparity of 0.00205 and a post-alignment
coordinate RMSE of 0.00193 across all 549 individuals, supporting robustness
of the leading two-dimensional population structure to high-LD-region masking.

## Pairwise FST

Genome-wide pairwise differentiation among CT, CEU, FIN, GBR, IBS and TSI is estimated with the Hudson FST estimator in PLINK2. The analysis uses the joint MAF>=0.05, high-LD-masked and LD-pruned PCA panel (287,815 SNPs). Standard errors are obtained by block jackknife using 500-SNP blocks.

FST is used as a descriptive complement to PCA. The former average-linkage clustering of pairwise FST values is not retained as a paper analysis.

## Runs of homozygosity

Primary ROH analysis:

- joint CT + EUR panel: 549 individuals
- joint MAF >=0.05: 4,878,327 SNPs
- no LD pruning
- PLINK `--homozyg`
- 50-SNP sliding windows
- minimum 50 SNPs per called ROH
- minimum length 1.5 Mb
- maximum density 50 kb/SNP (empirically non-binding; observed maximum 11.383 kb/SNP)
- maximum internal gap 250 kb
- <=1 missing call per 50-SNP window
- <=2 heterozygous calls per 50-SNP window
- window threshold 0.05
- no global `--homozyg-het` cap

Primary fROH denominator:

- 2.77e9 bp

Final production:

- 3,933 ROH >=1.5 Mb
- CT mean total ROH burden: 40.37 Mb
- CT mean FROH >=1.5 Mb: 0.01457
- CT mean ROH burden >=5 Mb: 18.53 Mb
- CT mean FROH >=5 Mb: 0.00669

Parameter choice was supported by explicit sensitivity analyses of missing
genotypes, heterozygote tolerance, maximum internal gap, MAF scope and LD
pruning. The >=5-Mb burden is retained as a secondary length-specific endpoint.

## LD decay

LD decay is estimated separately within each population using hard-call r2.

Primary decay definition:

`d_LD,p = min{d: mean_r2_p(d) < (1/3) * mean_r2_p(1 kb)}`

Operational details:

- 0.5-1.5 kb baseline
- 1-kb distance bins
- deterministic physical anchors
- first below-threshold bin = primary estimate
- five consecutive below-threshold bins = stability diagnostic only

For CT:

- MAF >=0.05 panel: 5,007,326 SNPs
- primary LD scale: 55.5 kb
- pre-phasing median: 99 SNPs

## Phasing

CT is phased with SHAPEIT2 v2.r904 using the 1000 Genomes Phase 3 EUR reference.

Settings:

- 503 EUR reference samples / 1006 haplotypes
- Ne=11,418
- window=0.5 Mb
- states=400
- burn/prune/main=7/8/20
- one thread per chromosome
- seed=15052011

Formal two-stage reference alignment:

- starting SNPs: 5,007,326
- absent from reference: 90,392
- allele-misaligned: 2,651
- retained: 4,914,283
- retained fraction: 98.14%

Post-phasing audit:

- 22 autosomes
- 46 samples
- 92 haplotypes
- dosage mismatches: 0
- duplicate sites: 0
- non-binary alleles: 0
- ordering failures: 0

## Selection

Primary statistic:

- saltiLASSI Lambda

Secondary statistic:

- original LASSI T

Implementation:

- `lassip` v1.2.1
- commit `a6a9d18c2323330fbf74d5a490f9e9c4ebe41d7c`
- phased input
- K=10
- Model D
- `--lassi-choice 4`
- physical distance
- `--max-extend-bp 100000`

Production populations:

- CT
- CEU
- TSI
- IBS

Production marker policy:

- autosomal biallelic SNPs
- within-population MAF >=0.05
- A/T and C/G excluded
- no LD pruning

CT exact phased-input density:

- physical scale: 55.5 kb
- median: 97 SNPs
- current saltiLASSI step candidate: 49 SNPs

CEU/TSI/IBS windows are recalibrated on the final non-palindromic production panels before scanning.

A separate CEU paper-oriented LASSI-T benchmark is kept distinct from the production CEU scan.

## Candidate regions

CT uses an empirical-outlier framework:

- primary threshold: top 1% genome-wide Lambda
- descriptive thresholds: top 0.1% and top 5%
- merge consecutive above-threshold windows

Candidate regions are interpreted as empirical outliers, not as simulation-based significant loci.

## CRG100

Production resource:

- UCSC hg19 `wgEncodeCrgMapabilityAlign100mer.bigWig`
- threshold: mean CRG100 >=0.9

Historical exploratory analyses used CRG >0.8 as the operational HIGHCONF cutoff and 0.9 as CRG_STRONG.

The production interval used for the mean CRG100 calculation is still being verified.
