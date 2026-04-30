# Parameter rationale

## LD pruning

Two different LD pruning thresholds are used:

- PCA: 50 5 0.2  
  Reduces correlation between SNPs and avoids distortion of principal components.

- ROH: 50 5 0.5  
  Preserves long homozygous segments while reducing local LD.

## MAF filtering

- MAF ≥ 0.05 used for:
  - PCA
  - Phasing
  - Selection analysis

This removes rare variants that may introduce noise and instability.

## HWE filtering

- Applied per chromosome  
- Bonferroni-corrected threshold  

Ensures removal of variants deviating from equilibrium due to technical artifacts.

## LD decay

Used to estimate the genomic scale of LD in the dataset.

This parameter is then used to define the window size in LASSI analyses.
