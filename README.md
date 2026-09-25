# Cinque Terre Population Genomics & Selection Pipeline

Reproducible Snakemake workflow for population-genomic analyses and haplotype-based selection scans in the Cinque Terre WGS cohort.

## Target dataset

The target dataset is the Cinque Terre WGS cohort.

- 50 WGS samples (all males)
- 11,490,646 SNPs
- Chromosomes `1-22` for downstream analyses
- GRCh37.p13 / hs37d5

The workflow starts from 50 individuals. Sample-level QC reproduces and documents the exclusion of four reviewed outliers, leaving **46 samples** for downstream population-genomic analyses.

Variant calling, joint genotyping, SNP/biallelic selection, and VQSR were performed upstream and are not repeated within this workflow. Their provenance is documented in `docs/upstream_variant_calling_provenance.md`.

## Target QC

The final autosomal CT core-QC dataset contains:

- 46 samples
- 9,000,246 SNPs
- strand-ambiguous A/T and C/G SNPs removed
- site missingness <=5%
- individual missingness <=5%
- chromosome-wise HWE filtering with Bonferroni threshold `0.05 / n_tested_chr`

No additional arbitrary QUAL, DP, GQ, or allele-balance hard filters are imposed after the reviewed sample exclusions.

## 1000 Genomes EUR reference

The 1000 Genomes Phase 3 EUR panel contains 503 individuals:

- CEU: 99
- FIN: 99
- GBR: 91
- IBS: 107
- TSI: 107

The panel is used for European context, CT phasing support, and population-specific selection comparisons.

For production selection scans:

- CEU, TSI, and IBS are comparison populations
- FIN and GBR are LD-context populations only
- production marker panels use within-population MAF >=0.05
- A/T and C/G SNPs are removed

CEU also has a separate paper-oriented LASSI-T benchmark, distinct from the production CEU comparison.

## Population-genomic analyses

Implemented analyses include:

- KING relatedness
- joint CT + EUR PCA
- runs of homozygosity and fROH
- ROH burden and length distribution
- population-specific LD decay
- SHAPEIT2 phasing
- saltiLASSI / LASSI selection scans

## Phasing

CT is phased with SHAPEIT2 v2.r904 using the 1000 Genomes Phase 3 EUR reference.

- CT MAF >=0.05 input: 5,007,326 SNPs
- SNPs retained after SHAPEIT2 reference alignment: 4,914,283
- retained samples: 46
- phased haplotypes: 92
- post-phasing genotype audit: PASS

## Selection scan

Primary method:

- **saltiLASSI Lambda**

Secondary method:

- **original LASSI T**

Implementation:

- `lassip` v1.2.1
- commit `a6a9d18c2323330fbf74d5a490f9e9c4ebe41d7c`
- K=10
- Model D / `--lassi-choice 4`
- physical-distance mode
- `--max-extend-bp 100000`

For CT, the primary LD scale is 55.5 kb. On the exact phased scan panel this corresponds to a median of **97 SNPs**, giving a current saltiLASSI step candidate of **49 SNPs**.

CEU, TSI, and IBS windows are being recalibrated after applying the same A/T and C/G exclusion used for CT.

The CT primary candidate threshold is the top 1% of genome-wide Lambda. Top 0.1% and top 5% are descriptive.

## Mappability

The production CRG100 resource is:

`wgEncodeCrgMapabilityAlign100mer.bigWig`

- build: hg19 / GRCh37
- expected MD5: `a1b1a8c99431fedf6a3b4baef028cca4`
- production threshold: mean CRG100 >=0.9

The exact interval used to calculate the mean CRG100 is being verified before the production filter is encoded.

## Repository structure

```text
config/       pipeline configuration
docs/         methods and provenance
workflow/     Snakemake rules and scripts
envs/         software environments
results/      generated outputs
logs/         execution logs
benchmarks/   runtime benchmarks
resources/    external project resources
```

Machine-specific paths belong in `config/config.local.yaml`, which is not version-controlled.

## Documentation

- `docs/upstream_variant_calling_provenance.md`
- `docs/target_qc_plan.md`
- `docs/methods_overview.md`
- `docs/pipeline_diagram.md`
- `docs/selection_strategy.md`
- `docs/selection_saltilassi_methodology.md`
- `docs/analysis_records/`

Historical exploratory runs are retained as provenance only and do not override the current workflow.
