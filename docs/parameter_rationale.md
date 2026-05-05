# Parameter rationale

## General principles

Parameter selection was guided by commonly adopted practices in population genomics and by the need to balance data quality, statistical power, and biological interpretability. Whenever possible, thresholds were chosen based on literature standards and adapted to the characteristics of the dataset.


---

## Variant-level QC

### Missingness (geno, mind)

A threshold of 0.05 was applied for both variant-level (geno) and sample-level (mind) missingness. This value represents a commonly used compromise between retaining sufficient data and excluding poorly genotyped variants or samples that could introduce bias (Purcell et al., 2007).

Stricter thresholds (e.g., 0.03) may be explored in sensitivity analyses.


### Hardy–Weinberg equilibrium

HWE filtering was performed using a Bonferroni-corrected threshold (α = 0.05), applied separately for each chromosome. This approach controls for multiple testing and reduces the inclusion of variants affected by genotyping errors or strong deviations from equilibrium (Purcell et al., 2007).


---

## Relatedness

Two complementary approaches were used:

- KING, which does not require LD pruning and is robust for detecting close relationships  
- PLINK IBS (--genome), applied on LD-pruned common variants  

The use of both methods allows cross-validation of relatedness estimates and increases reliability in identifying duplicates and relatives.


---

## PCA

LD pruning parameters (50 SNP window, 5 SNP step, r² = 0.2) were selected to remove correlated variants and reduce redundancy while preserving genome-wide structure.

A MAF threshold of 0.05 was applied to focus on common variants that capture population-level structure more effectively (Price et al., 2006).


---

## ROH

A less stringent LD pruning threshold (r² = 0.5) was used compared to PCA. This choice preserves more local genomic structure, which is important for detecting runs of homozygosity.

ROH segments were classified based on their length distribution to distinguish between recent and ancient inbreeding signals (Ceballos et al., 2018).


---

## LD decay

LD decay was estimated to determine the genomic scale of correlation between variants. The distance at which LD decays to one-third of its initial value was used to define the window size for LASSI analyses.

This approach allows data-driven parameterization rather than relying on arbitrary window sizes.


---

## Phasing

Only SNPs with MAF ≥ 0.05 were retained for phasing to improve accuracy, as rare variants are more prone to phasing errors and contribute less to haplotype-based inference.

Reference populations (CEU, TSI, IBS) were selected to match the expected ancestry of the target population (Delaneau et al., 2012).


---

## LASSI

The likelihood window size and step were derived from LD decay estimates. This ensures that the window captures meaningful haplotype structure while avoiding excessive smoothing.

The shift was set to 10% of the window size, following common practice in sliding-window analyses.


---

## Interpretation of selection signals

The LASSI m parameter was used to distinguish between:

- Hard sweeps (m = 1)  
- Soft sweeps (m > 1)  

This classification provides insight into the mode of selection acting on candidate regions (Abondio et al., 2022).


---

## Optional demographic analysis

IBDNe was included as an optional analysis to infer recent effective population size (Ne).

This method relies on the distribution of IBD segment lengths, where longer segments reflect recent shared ancestry and shorter segments capture older demographic events.


---

## References

Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A. R., Bender, D., Maller, J., Sklar, P., de Bakker, P. I. W., Daly, M. J., & Sham, P. C. (2007). PLINK: A tool set for whole-genome association and population-based linkage analyses. *The American Journal of Human Genetics, 81*(3), 559–575. https://doi.org/10.1086/519795

Price, A. L., Patterson, N. J., Plenge, R. M., Weinblatt, M. E., Shadick, N. A., & Reich, D. (2006). Principal components analysis corrects for stratification in genome-wide association studies. *Nature Genetics, 38*(8), 904–909. https://doi.org/10.1038/ng1847

Ceballos, F. C., Joshi, P. K., Clark, D. W., Ramsay, M., & Wilson, J. F. (2018). Runs of homozygosity: Windows into population history and trait architecture. *Nature Reviews Genetics, 19*(4), 220–234. https://doi.org/10.1038/nrg.2017.109

Delaneau, O., Marchini, J., & Zagury, J. F. (2012). A linear complexity phasing method for thousands of genomes. *Nature Methods, 9*(2), 179–181. https://doi.org/10.1038/nmeth.1785
