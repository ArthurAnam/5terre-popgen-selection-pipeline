# 5Terre Population Genomics & Selection Pipeline

Reproducible workflow for population structure and natural selection analyses in the Cinque Terre population.

---

## 🧬 Overview

This project implements a modular and reproducible pipeline for population genomics analyses, including:

- Variant-level quality control (QC)
- Relatedness filtering (IBS / KING)
- Population structure (PCA)
- Runs of Homozygosity (ROH)
- Linkage Disequilibrium (LD) decay
- Haplotype phasing (SHAPEIT2)
- Selection scans (LASSI)
- Functional enrichment analysis

All analyses are designed to be fully reproducible using Snakemake.

---

## 🧠 Study design

The pipeline is built around a branch-based design:

- **QC backbone** → common preprocessing  
- **Population structure branch** → PCA, ROH  
- **Selection branch** → LD decay → phasing → LASSI → interpretation  
- **Optional demography branch** → IBD / IBDNe  

Each analysis operates on a dataset specifically tailored for its purpose.

---

## ⚙️ Technologies

- Snakemake (workflow manager)
- PLINK
- bcftools
- SHAPEIT2
- EIGENSOFT (smartpca)
- LASSI

---

## 📁 Repository structure


config/ → pipeline configuration
docs/ → methods, rationale, diagrams
workflow/ → Snakemake rules (in development)
envs/ → environments (conda)
results/ → outputs (generated)
logs/ → logs
benchmarks/ → runtime benchmarks


---

## 📄 Documentation

- Pipeline design: `docs/pipeline_diagram.md`
- Methods overview: `docs/methods_overview.md`
- Parameter rationale: `docs/parameter_rationale.md`
- Reproducibility principles: `docs/reproducibility_principles.md`

---

## 🚧 Status

This repository currently focuses on pipeline design and reproducibility structure.  
Execution rules and datasets will be integrated in future iterations.

---

## 📌 Notes

- Genome build: GRCh37.p13 (hg19)
- SNP-only, biallelic, autosomal variants
- Branch-specific preprocessing strategy

---

## 👤 Author

Tommaso  Barbaresi
