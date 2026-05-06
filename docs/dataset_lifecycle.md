# Dataset lifecycle

## Overview

This document describes the lifecycle of datasets throughout the pipeline, including their origin, preprocessing steps, intended analytical use, and relationships between analysis branches.

The pipeline distinguishes between:

- raw datasets
- quality-controlled datasets
- harmonized datasets
- branch-specific datasets

This separation is intended to preserve reproducibility, traceability, and biological interpretability.

---

## Original datasets

### Cinque Terre

Target population under study.

Expected data type:

- Whole-genome sequencing (WGS)
- Joint-called VCF

Role in the pipeline:

- Main target population
- Used for population structure analyses
- Used for selection analyses
- Used for optional demographic inference

Expected characteristics:

- Human diploid genomes
- GRCh37.p13 / hg19 coordinates
- SNP-focused analyses
- Autosomal, biallelic SNPs retained for downstream analyses

---

### 1000 Genomes EUR

Reference dataset used for:

- European ancestry contextualization
- PCA comparison
- Reference-guided phasing

Current populations:

- CEU
- TSI
- IBS

Additional EUR populations potentially available:

- FIN
- GBR

Expected characteristics:

- Already phased
- GRCh37 / hg19 coordinates
- External reference-only role

The 1000 Genomes dataset is not intended to be directly merged into the selection scan branch after phasing.

---

### Ferrara cohorts

Future external cohorts potentially including:

- Northern Italy
- Corsica
- Sardinia
- France

Expected role:

- Comparative population structure analyses
- Contextualization of the Cinque Terre population

Current status:

- Not yet integrated

Expected strategy:

- Harmonized with the main dataset only within the population structure branch
- Not directly incorporated into the LASSI selection branch

---

## Dataset preprocessing lifecycle

### RAW_DATASETS

Initial unmodified datasets.

Characteristics:

- Original files
- No filtering
- No harmonization
- No normalization

Storage location:

- `data/raw/`

Purpose:

- Permanent source datasets
- Reproducibility reference
- Never modified directly

---

### HARMONIZED_DATASETS

Datasets generated after technical harmonization.

Potential harmonization steps:

- chromosome naming consistency
- coordinate/build check
- reference allele consistency check
- duplicate variant removal
- retention of autosomal variants
- retention of biallelic SNPs only

Purpose:

- Generate technically comparable datasets before downstream QC and merging

Storage location:

- `data/processed/`

---

### QC_CLEAN

Dataset generated after variant-level quality control.

Main filtering strategy:

- Autosomal variants only
- SNP-only dataset
- Biallelic variants only
- Multiallelic variants excluded
- Monomorphic variants removed
- Palindromic SNPs removed
- Missingness filtering
- HWE filtering

Purpose:

- Shared high-quality dataset for downstream analyses

Used by:

- Relatedness filtering
- Population structure branch
- Selection branch

---

### MASTER_UNRELATED

Dataset generated after relatedness filtering.

Derived from:

- QC_CLEAN

Relatedness analyses:

- KING
- PLINK IBS / `--genome`

Purpose:

- Reduce confounding from close relatives
- Generate a dataset suitable for population structure analyses

Used by:

- PCA
- ROH
- Harmonized population structure analyses

---

## Branch separation strategy

The pipeline separates population structure analyses from selection analyses after initial quality control.

This separation is intended to preserve biologically meaningful haplotype structure while minimizing confounding effects introduced by dataset merging.

---

## Population structure branch

Datasets potentially included:

- Cinque Terre
- 1000 Genomes EUR
- Future Ferrara cohorts

Expected preprocessing:

- Dataset harmonization
- Shared SNP intersection
- LD pruning
- Potential merge operations

Main analyses:

- PCA
- ROH
- Relatedness exploration
- Population contextualization

Purpose:

- Characterize ancestry and population structure
- Compare the target population with external references

---

## Selection branch

Primary dataset:

- Cinque Terre

Reference support:

- 1000 Genomes EUR phased haplotypes

Main preprocessing:

- MAF filtering
- LD decay estimation
- Phasing

Main analyses:

- LASSI
- Sweep classification
- Candidate region annotation
- Functional enrichment

Important strategy:

The target population is intentionally kept separate from external cohorts during selection scans in order to preserve population-specific haplotype structure.

---

## Optional demographic branch

Optional analyses based on:

- Phased haplotypes
- IBD segment detection

Potential tools:

- hap-IBD
- Refined IBD
- IBDNe

Purpose:

- Infer recent effective population size (Ne)
- Characterize recent demographic history

Expected outputs:

- Ne trajectories over time
- IBD segment distributions

Demographic estimates are interpreted in generations and may be converted into calendar years using a generation time of 29 years.

---

## Final outputs

Expected final outputs include:

- PCA coordinates and figures
- ROH metrics and distributions
- LD decay curves
- Phased haplotypes
- LASSI likelihood profiles
- Candidate sweep regions
- Functional enrichment analyses
- Optional demographic inference outputs

Output location:

- `results/`
