# Methods overview

## Study design

This study investigates population structure and signatures of natural selection in the 5Terre population using genome-wide SNP data. The analysis is based on a reproducible workflow implemented with Snakemake, ensuring transparency, modularity, and scalability.

All analyses are performed on SNP-only, biallelic, autosomal variants mapped to the GRCh37.p13 (hg19) reference genome.


---

## Variant-level quality control

A common variant-level quality control (QC) procedure is applied to all datasets prior to downstream analyses. The following filters are implemented:

- Retention of autosomal variants only  
- Inclusion of SNPs only  
- Removal of multiallelic variants  
- Removal of monomorphic variants  
- Removal of palindromic SNPs (A/T and C/G)  

Missingness filters are applied at both variant and sample levels:

- Variant missingness threshold: geno ≤ 0.05  
- Sample missingness threshold: mind ≤ 0.05  

Hardy–Weinberg equilibrium (HWE) filtering is performed separately for each chromosome using a Bonferroni-corrected significance threshold (α = 0.05).


---

## Relatedness and sample filtering

Sample-level quality control is performed to identify duplicates and related individuals.

Two complementary approaches are used:

- KING, applied without LD pruning and without MAF filtering  
- PLINK IBS estimation (--genome), applied on LD-pruned variants (window = 50 SNPs, step = 5 SNPs, r² threshold = 0.2) with MAF ≥ 0.05  

Based on these analyses, a filtered dataset of unrelated individuals (MASTER_UNRELATED) is defined and used for downstream population structure analyses.


---

## Population structure analyses

Population structure is investigated using the MASTER_UNRELATED dataset.

### Principal Component Analysis (PCA)

PCA is performed using smartpca (EIGENSOFT) after:

- LD pruning (50 SNP window, 5 SNP step, r² = 0.2)  
- MAF filtering (MAF ≥ 0.05)  

The analysis aims to:

- Assess ancestry of the 5Terre population  
- Position samples within the European genetic landscape  
- Compare the target population with reference populations from the 1000 Genomes Project (EUR panel)  


### Runs of Homozygosity (ROH)

ROH are computed using PLINK to characterize autozygosity patterns:

- SNP-only dataset  
- MAF ≥ 0.05  
- LD pruning optimized for ROH detection (50 SNP window, 5 SNP step, r² = 0.5)  

ROH metrics are used to infer:

- Individual autozygosity burden  
- Distribution of ROH segment lengths  
- Demographic history signals (e.g., isolation, inbreeding)  


---

## Linkage disequilibrium decay

LD decay is estimated using combined datasets (5Terre + 1000 Genomes EUR) to determine the genomic scale of correlation between variants.

The decay of LD is evaluated as a function of physical distance, and the characteristic distance at which LD reaches one-third of its initial value is estimated.

This value is used to define the window size for downstream selection analyses.


---

## Phasing

Haplotype phasing is performed on the 5Terre dataset using SHAPEIT2.

Reference haplotypes from selected European populations (CEU, TSI, IBS) from the 1000 Genomes Project are used as a reference panel.

Only SNPs with MAF ≥ 0.05 are retained for phasing.


---

## Detection of selection signals (LASSI)

Selection scans are performed using LASSI, a likelihood-based method for detecting selective sweeps from haplotype data.

Input data consist of phased haplotypes from the 5Terre population.

The likelihood window size is defined based on LD decay estimates, and a sliding window approach is used across the genome.

The method identifies candidate regions under selection and assigns a likelihood-based score to each genomic window.


---

## Post-processing and annotation

Candidate regions identified by LASSI are further processed and annotated:

- Aggregation of overlapping windows  
- Ranking of signals based on likelihood scores  
- Classification of regions as genic or intergenic  
- Annotation of overlapping genes  

Selective sweeps are classified based on the LASSI m parameter:

- Hard sweeps: m = 1  
- Soft sweeps: m > 1  


---

## Functional enrichment analysis

Genes overlapping candidate regions are analyzed using:

- g:Profiler (functional enrichment)  
- STRING (protein–protein interaction networks)  

These analyses provide biological interpretation of the detected selection signals.


---

## Optional demographic analysis

An optional analysis of recent demographic history can be performed using IBDNe.

This requires:

- Phased genotype data  
- Detection of IBD segments between individuals  

IBDNe estimates effective population size (Ne) over time based on the distribution of IBD segment lengths.
