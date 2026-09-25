# Cinque Terre Population Genomics & Selection Pipeline

Reproducible workflow for population structure analyses and haplotype-based selection scans focused on the Cinque Terre population.

---

## 🧬 Overview

This repository contains a modular and reproducible population genomics pipeline designed to support:

* input verification and population-genomics QC
* dataset harmonization
* population structure analyses
* relatedness and IBS exploration
* runs of homozygosity analyses
* LD decay estimation
* haplotype phasing with SHAPEIT2
* LASSI-based selection scans
* candidate region definition and prioritization
* functional genomic annotation and exploratory biological contextualization

The workflow is designed for reproducible execution using Snakemake and centralized configuration files.

---

## 🧠 Study design

The pipeline is organized around four main components:

* **Input verification and population-genomics QC**  
  Validate the delivered Cinque Terre WGS VCF, then apply downstream sample/variant QC without repeating upstream variant calling or VQSR.

* **Dataset harmonization**  
  Harmonize the QCed Cinque Terre target dataset with Italian WGS reference cohorts and 1000 Genomes Project EUR reference data where required.

* **Population structure analyses**  
  PCA, runs of homozygosity, fROH, ROH burden, ROH length distribution, and IBS / pairwise relatedness exploration.

* **Selection scan and candidate region interpretation**  
  LD decay estimation, phasing, LASSI scan, candidate region definition, ranking, and functional genomic contextualization.

An optional demographic branch based on IBD segment detection and IBDNe may be added depending on data suitability and final analytical priorities.

---

## 📊 Target dataset

The target dataset is the Cinque Terre WGS cohort.

**Verified starting dataset (`all_samples_snp.vcf.gz`):**

* 50 WGS samples
* 11,490,646 SNP records
* 0 INDELs
* 0 multiallelic sites
* all records `FILTER=PASS`
* indexed contigs `1-22`, `X`, `Y`
* GRCh37.p13 / hs37d5 provenance upstream

The workflow starts from all 50 delivered individuals. Sample-level QC then reproduces and documents the exclusion of four reviewed outliers, leaving **46 samples** for downstream population-genomic analyses.

Variant calling, joint genotyping, SNP/biallelic selection, and VQSR were performed upstream and are not repeated within this workflow. Their provenance is documented in `docs/upstream_variant_calling_provenance.md`.

---

## 🌍 Reference and comparison datasets

The main comparative framework is expected to involve Italian WGS reference cohorts, ideally representing:

* North Italy
* South Italy

These cohorts are currently under evaluation and may be based on published Italian WGS data, including Sazzini et al. (2020), depending on data availability and suitability.

The 1000 Genomes Project EUR panel is used as European reference context and phasing support.

Candidate 1000 Genomes EUR populations currently under consideration for selection-related comparison include:

* TSI
* FIN
* IBS

The final comparison set remains under evaluation.

---

## ⚙️ Technologies

Main tools and frameworks include:

* Snakemake
* PLINK
* KING
* bcftools
* SHAPEIT2
* EIGENSOFT / smartpca
* LASSI
* hap-IBD / Refined IBD, if optional demographic analyses are performed
* IBDNe, if optional demographic analyses are performed

The main Conda environment is defined in `envs/pipeline.yaml`.

---

## 📁 Repository structure

```text
config/         → central pipeline configuration
docs/           → methods, provenance, workflow design, and software framework
workflow/       → Snakemake rules
envs/           → Conda environment definition
data/           → input and processed data placeholders
results/        → generated outputs
logs/           → execution logs
benchmarks/     → runtime benchmarks
```

Exploratory QC rules remain available for targeted reruns, but their outputs are not part of the default `rule all`. The production target reproduces the final analysis from the reviewed, version-controlled QC decisions.

Generated outputs are not intended to be manually edited.

---

## 📄 Documentation

Current documentation:

* Upstream variant-calling provenance: `docs/upstream_variant_calling_provenance.md`
* Methods overview: `docs/methods_overview.md`
* Pipeline diagram: `docs/pipeline_diagram.md`
* Selection strategy: `docs/selection_strategy.md`
* Computational framework: `docs/computational_framework.md`

Central configuration:

* `config/config.yaml`

Machine-specific input paths:

* `config/config.local.yaml` (local only; intentionally not version-controlled)

Main Conda environment:

* `envs/pipeline.yaml`

---

## 🚧 Status

The reproducible pre-QC stage is implemented and has validated the current 50-sample PASS-only SNP input dataset. Downstream QC and analysis modules will be integrated incrementally and checked at each stage.

---

## 📌 Notes

Current working assumptions and principles:

* Genome build: GRCh37.p13 / hg19 coordinates, with hs37d5 upstream provenance
* Coordinate system: VCF 1-based
* Downstream population-genomics analyses: autosomal chromosomes 1-22 unless otherwise specified
* Delivered variant class: biallelic PASS SNPs, independently verified by pre-QC
* Historical 46-sample result: to be reproduced rather than assumed
* HWE filtering strategy is under methodological review before being applied to this small target population
* Selection results are interpreted as candidate genomic regions deviating from neutral expectations, not as definitive proof of adaptive loci or causal variants
* Functional interpretation is treated as exploratory biological contextualization
