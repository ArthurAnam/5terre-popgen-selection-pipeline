# Target WGS QC plan

The delivered Cinque Terre VCF contains 50 WGS samples and has already undergone upstream variant calling and VQSR. Four samples have been reviewed and approved for exclusion after sample-level QC: TSBC6060, TSBC6389, TSBD8047 and TSBD8199. Variant-level QC will start from the reviewed 46-sample dataset.

## Principle

Diagnostic exploration and hard filtering are kept distinct. Classical site/genotype quality annotations are first described and plotted. They do not automatically become additional hard filters because the input has already passed GATK VQSR, which uses a multivariate quality model. A new hard threshold is added only if the post-VQSR data show a specific residual problem that justifies it.

## Variant-level quantities explored in this pipeline

QUAL is inspected descriptively. MAF and genotype missingness are recalculated on the retained post-sample-QC cohort. Variant category counts include monomorphic sites, palindromic A/T and C/G SNPs, multiallelic/non-SNP records, and sites above the planned 5% missingness threshold.

VQSR-specific bookkeeping fields (VQSLOD, culprit, POSITIVE_TRAIN_SITE and NEGATIVE_TRAIN_SITE) are not used for downstream QC because VQSR was completed upstream. The pipeline does not attempt to reproduce or second-guess the upstream VQSR model.

Upstream annotations such as QD, FS, SOR, MQ and rank-sum metrics remain available in the delivered VCF but are not turned into new hard filters by default.

## Available genotype-level annotations to explore

FORMAT/DP, FORMAT/GQ, FORMAT/AD and GT.

These will be summarized descriptively before any decision about genotype-level masking or hard cutoffs.

## Core filtering decisions

1. Remove the four reviewed sample-level QC failures.
2. Restrict to autosomes 1-22.
3. Retain SNPs and biallelic variants only. This is expected to be a no-op for the delivered VCF but remains an explicit validation step.
4. Remove variants that become monomorphic after sample removal.
5. Remove strand-ambiguous / palindromic SNPs (A/T and C/G) as part of the core target QC.
6. Remove variants with missingness >5%.
7. Remove individuals with missingness >5%.
8. Recompute site summaries after any sample removal and remove newly monomorphic variants if necessary.
9. Recompute HWE on the retained sample set and apply the agreed HWE rule.

## Exploration before any additional hard filter

QUAL and current cohort-level missingness/frequency summaries are explored first without filtering. FORMAT-level DP, GQ and AD are explored separately because they describe genotype calls. No VQSR-specific field is used as a new downstream filter, and no external threshold is copied automatically.

## HWE

HWE is part of the intended QC and must be recomputed after sample-level QC rather than inferred from upstream INFO annotations. The current configuration specifies a chromosome-wise Bonferroni threshold. Before applying the hard filter, report how many sites would be removed and their distribution so that the impact of the choice is transparent.

## Missingness

The planned threshold is 5% for both variants and individuals. Missingness should be evaluated iteratively because sample removal can alter site missingness and create new monomorphic variants.

## MAF

No general MAF filter is part of the core WGS QC. MAF thresholds remain analysis-specific for later branches.
