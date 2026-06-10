# Cinque Terre Population Genomics & Selection Pipeline

Reproducible workflow for population structure analyses and haplotype-based selection scans in the Cinque Terre population.

---

## 🧬 Overview

This repository contains a modular and reproducible population genomics pipeline designed to support:

* dataset harmonization and population-genomics QC
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

* **Input datasets and harmonization**
  Cinque Terre WGS data, Italian WGS reference cohorts under evaluation, and 1000 Genomes Project EUR reference data.

* **Population-genomics QC**
  Autosomal biallelic SNP filtering, missingness filtering, strand-ambiguous SNP removal, and chromosome-wise HWE filtering using Bonferroni correction.

* **Population structure analyses**
  PCA, runs of homozygosity, fROH, ROH burden, ROH length distribution, and IBS / pairwise relatedness exploration.

* **Selection scan and candidate region interpretation**
  LD decay estimation, phasing, LASSI scan, candidate region definition, ranking, and functional genomic contextualization.

An optional demographic branch based on IBD segment detection and IBDNe may be added depending on data suitability and final analytical priorities.

---

## 📊 Target dataset

The target dataset is the Cinque Terre WGS cohort.

Expected characteristics:

* 46 WGS samples
* mean depth approximately 30×
* GRCh37.p13 / hg19 coordinates
* variants already called upstream
* variants already VQSR-filtered upstream

Variant calling and VQSR are not repeated within this workflow.

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

---

## 📁 Repository structure

```text
config/         → central pipeline configuration
docs/           → methods, workflow design, selection strategy, and software framework
workflow/       → Snakemake rules
envs/           → Conda environments
data/           → input and processed data placeholders
results/        → generated outputs
logs/           → execution logs
benchmarks/     → runtime benchmarks
```

Generated outputs are not intended to be manually edited.

---

## 📄 Documentation

Current documentation:

* Methods overview: `docs/methods_overview.md`
* Pipeline diagram: `docs/pipeline_diagram.md`
* Selection strategy: `docs/selection_strategy.md`
* Computational framework: `docs/computational_framework.md`

Central configuration:

* `config/config.yaml`

---

## 🚧 Status

This repository currently focuses on pipeline design, documentation, and reproducibility structure.

Execution rules, finalized datasets, and analysis-specific scripts will be integrated progressively.

---

## 📌 Notes

Current working assumptions:

* Genome build: GRCh37.p13 / hg19
* Coordinate system: VCF 1-based
* Main variant class: autosomal biallelic SNPs
* HWE filtering: chromosome-wise Bonferroni correction
* Selection results are interpreted as candidate genomic regions deviating from neutral expectations, not as definitive proof of adaptive loci or causal variants
* Functional interpretation is treated as exploratory biological contextualization
