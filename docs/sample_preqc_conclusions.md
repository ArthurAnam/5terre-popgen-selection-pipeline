# Individual sample pre-QC conclusions

## Scope

This document records the reviewed sample-level QC decision for the 50-sample
Cinque Terre WGS input dataset. The decision was made from the current
post-VQSR VCF and was not based on a target post-QC sample count.

All diagnostic comparisons were restricted to autosomes 1-22. Statistical
flags were used to identify samples for review; no sample was excluded
automatically from a single threshold.

## Evidence reviewed

The initial autosomal sample pre-QC compared average depth, genotype
missingness/call rate and heterozygosity as primary technical metrics, with
non-reference variant count, singleton fraction and per-sample Ti/Tv as
secondary diagnostics. Robust z-scores were calculated relative to the
50-sample cohort.

Chromosome-level follow-up showed that the principal abnormalities were
genome-wide rather than driven by isolated loci or chromosomes.

A second diagnostic step inspected VCF genotype fields GT, DP, GQ and AD over
a deterministic genome-wide autosomal subset. These metrics were compared
between samples without applying external DP, GQ or allele-balance cutoffs.

## Reviewed exclusions

### TSBC6060

- Excess heterozygosity was present on 22/22 autosomes.
- Average depth and missingness were not broadly abnormal.
- Heterozygous allele balance was strongly shifted relative to the cohort.
- The combination supports a sample-level technical artefact affecting
  genotype composition.

### TSBC6389

- Excess heterozygosity was present on 22/22 autosomes.
- Reduced depth and increased missingness were present on a subset of
  chromosomes.
- Mean genotype quality was lower than the cohort.
- Heterozygous allele balance was strongly shifted relative to the cohort.
- The combined evidence supports a sample-level technical artefact.

### TSBD8047

- Excess heterozygosity was present on 22/22 autosomes.
- Average depth and missingness were not broadly abnormal.
- Heterozygous allele balance was strongly shifted relative to the cohort.
- The combination supports a sample-level technical artefact affecting
  genotype composition.

### TSBD8199

- Low average depth was present on 22/22 autosomes.
- Elevated missingness was present on 22/22 autosomes.
- Heterozygosity was abnormal on 21/22 autosomes.
- Genotype-field diagnostics confirmed very low depth and substantially lower
  genotype quality relative to the cohort.
- This is a clear low-coverage / low-quality sample.

## Interpretation and terminology

TSBC6060, TSBC6389 and TSBD8047 are described as showing a pattern consistent
with a sample-level technical artefact. The present VCF-based QC does not
identify the exact cause. In particular, contamination or sample mixture
should not be stated as a confirmed diagnosis without a dedicated
contamination analysis, ideally using the underlying BAM/CRAM files.

TSBD8199 is classified as a low-coverage / low-quality sample.

The four samples are therefore approved for exclusion before variant-level QC:

- TSBC6060
- TSBC6389
- TSBD8047
- TSBD8199

This exclusion decision is explicit and reviewed. It is not implemented as an
automatic rule such as "remove every sample beyond a fixed robust-z
threshold".
