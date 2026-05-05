# Parameter rationale

## General principles

Parameter selection was guided by commonly adopted practices in population genomics and by the need to balance data quality, statistical power, and biological interpretability. Thresholds were selected to maximize robustness while preserving sufficient variation for downstream analyses.


---

## Variant-level QC

### Missingness (geno, mind)

A threshold of 0.05 was applied for both variant-level (geno) and sample-level (mind) missingness. This represents a balance between removing poorly genotyped data and retaining sufficient variants and samples.

Stricter thresholds may be explored in sensitivity analyses if needed.


### Hardy–Weinberg equilibrium

HWE filtering was performed using a Bonferroni-corrected threshold (α = 0.05), applied separately for each chromosome. This reduces the inclusion of variants affected by genotyping errors or strong deviations from equilibrium.


---

## Relatedness

Two complementary approaches were used:

- KING, which does not require LD pruning and is robust for detecting close relationships  
- PLINK IBS (--genome), applied on LD-pruned common variants  

Using both methods allows cross-validation of relatedness estimates and improves reliability.


---

## PCA

LD pruning parameters (50 SNP window, 5 SNP step, r² = 0.2) were selected to remove correlated variants while preserving genome-wide structure.

A MAF threshold of 0.05 was used to focus on common variants that better capture population-level structure.


---

## ROH

A less stringent LD pruning threshold (r² = 0.5) was used compared to PCA to preserve local genomic structure.

ROH segments will be classified based on their length distribution to distinguish between recent and ancient inbreeding signals.


---

## LD decay

LD decay is estimated to determine the genomic scale of correlation between variants.

The distance at which LD decays to one-third of its initial value is used to define the window size for LASSI analyses, enabling data-driven parameter selection.


---

## Phasing

Only SNPs with MAF ≥ 0.05 were retained for phasing to improve accuracy.

Reference populations (CEU, TSI, IBS) were selected to match the expected ancestry of the target population.


---

## LASSI

The likelihood window size and step are derived from LD decay estimates to capture meaningful haplotype structure.

The shift is set to 10% of the window size, following standard sliding-window approaches.


---

## Interpretation of selection signals

The LASSI m parameter is used to distinguish between:

- Hard sweeps (m = 1)  
- Soft sweeps (m > 1)  

This classification provides insight into the mode of selection acting on candidate regions.


---

## Optional demographic analysis

IBDNe is included as an optional analysis to infer recent effective population size (Ne).

This method relies on the distribution of IBD segment lengths, where longer segments reflect recent shared ancestry and shorter segments capture older demographic events.
