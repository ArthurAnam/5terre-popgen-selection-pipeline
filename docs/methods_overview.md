# Methods overview

The workflow is divided into three main analytical components:

## 1. Variant-level QC

All datasets are harmonized and filtered prior to downstream analyses:
- autosomal SNPs only
- removal of multiallelic variants
- removal of palindromic SNPs
- HWE filtering (Bonferroni-corrected threshold)

## 2. Population structure

Performed on LD-pruned datasets:
- PCA (smartpca)
- IBS / relatedness (PLINK / KING)
- ROH analysis

Different datasets are used for different analyses to avoid bias.

## 3. Selection analysis

Performed on non-pruned datasets:
- LD decay estimation
- MAF filtering (≥ 0.05)
- Phasing (SHAPEIT2)
- LASSI scan
- Annotation and enrichment

This ensures that haplotype-based methods are not affected by LD pruning.
