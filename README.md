# Cinque Terre Population Genomics & Selection Pipeline

Reproducible Snakemake workflow for population-genomic analyses and haplotype-based selection scans focused on the Cinque Terre WGS cohort.

---

## Overview

The repository contains the version-controlled analysis workflow for:

* input verification and population-genomics QC
* 1000 Genomes EUR harmonization and European reference analyses
* PCA, relatedness and runs of homozygosity
* LD-decay calibration
* reference-assisted SHAPEIT2 phasing of Cinque Terre
* saltiLASSI selection scans, with original LASSI T retained as a secondary benchmark
* candidate/outlier-region definition and downstream functional contextualization

Variant calling, joint genotyping, SNP/biallelic selection and VQSR were performed upstream and are not repeated here.

---

## Target dataset

**Verified delivered dataset (`all_samples_snp.vcf.gz`):**

* 50 WGS samples
* 11,490,646 SNPs
* 0 INDELs
* 0 multiallelic sites
* all records `FILTER=PASS`
* indexed contigs `1-22`, `X`, `Y`
* GRCh37.p13 / hs37d5 upstream provenance

Sample-level QC reproduces and documents the exclusion of four reviewed outliers, leaving **46 samples** for downstream analyses.

The final autosomal Cinque Terre core-QC dataset contains **9,000,246 SNPs**. The MAF>=0.05 selection/phasing input contains **5,007,326 SNPs** before the formal SHAPEIT2 reference-alignment check.

---

## Reference and comparison datasets

The current European reference framework uses the **503-sample 1000 Genomes Phase 3 EUR panel**:

* CEU: 99
* FIN: 99
* GBR: 91
* IBS: 107
* TSI: 107

For the production selection comparison, the frozen population set is **CEU, TSI and IBS**. FIN and GBR are retained as European LD context only.

Production selection panels are population-specific and use MAF>=0.05 with strand-ambiguous A/T and C/G SNPs excluded. CEU will additionally have a separate methodological benchmark designed to approximate the Harris & DeGiorgio (2020) empirical LASSI-T setup; that benchmark is not the production CEU comparator.

Italian WGS reference cohorts remain a possible future extension for broader population-structure comparisons and are not required for the current selection branch.

---

## Current analytical state

Completed and audited:

* target sample and variant QC
* chromosome-wise HWE filtering
* KING relatedness audit
* CT + 1000G EUR harmonization
* joint PCA and high-LD-region sensitivity
* ROH/fROH analyses
* CT and EUR LD-decay audits
* SHAPEIT2 production phasing of all 22 CT autosomes
* exhaustive post-phasing structural/genotype audit
* pre-LASSI environment and lassip software gates
* CRG100 resource gate
* exact CT phased-input audit for lassip

The CT phased lassip input contains **4,914,283 SNPs** across 46 individuals. Recounting SNP density at the primary 55.5-kb LD scale gives a median of **97 SNPs** (49-SNP half-window step candidate), so the older pre-phasing 99-SNP value is retained only as provenance until the production window is formally re-frozen.

CEU/TSI/IBS LD-derived SNP windows are being recalibrated after applying the same production A/T and C/G exclusion used for CT. No production saltiLASSI genome-wide scan has been started yet.

---

## Core QC principles

For Cinque Terre:

* autosomes 1-22
* A/T and C/G SNPs removed
* monomorphic sites removed
* site missingness >5% removed
* individual missingness >5% evaluated iteratively with cohort-dependent site metrics recomputed after sample removal
* chromosome-wise HWE Bonferroni threshold: `0.05 / n_tested_chr`
* no additional arbitrary QUAL, DP, GQ or allele-balance hard filters after the reviewed sample exclusions

The delivered VCF was already restricted upstream to PASS biallelic SNPs; this is verified, not claimed as a new downstream filtering step.

---

## Selection design

Primary method:

* **saltiLASSI Lambda** using `lassip` v1.2.1

Secondary method:

* **original LASSI T** from the same maintained implementation

Frozen model-level settings include phased input, K=10, Model D (`--lassi-choice 4` in lassip v1.2.1), physical-distance mode and `--max-extend-bp 100000`.

CT inference uses an empirical-outlier framework: top 1% genome-wide Lambda is the primary candidate threshold; top 0.1% and top 5% are descriptive. Consecutive above-threshold windows are merged into candidate/outlier regions. These are not described as simulation-based significant regions.

The production CRG100 threshold is 0.9. The exact interval/unit on which mean CRG100 is applied is being verified before the production filter is encoded.

---

## Repository structure

```text
config/         central configuration
docs/           methods, provenance and analysis records
workflow/       Snakemake rules and scripts
envs/           software environment definitions
results/        generated outputs
logs/           execution logs
benchmarks/     runtime/resource benchmarks
resources/      project-contained external resources where appropriate
```

Machine-specific paths belong in `config/config.local.yaml`, which is intentionally not version-controlled. Large genomic resources and generated results are not committed.

---

## Main documentation

* `docs/upstream_variant_calling_provenance.md`
* `docs/target_qc_plan.md`
* `docs/methods_overview.md`
* `docs/pipeline_diagram.md`
* `docs/selection_strategy.md`
* `docs/selection_saltilassi_methodology.md`
* `docs/pre_lassi_environment_audit.md`
* `docs/analysis_records/`

Historical exploratory outputs are retained only as provenance and must not silently override the current post-16-September-2026 workflow decisions.
