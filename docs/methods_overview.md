# Methods

## Study design

This study investigates population structure and signatures of natural selection in the 5Terre population using genome-wide SNP data. All analyses are performed on SNP-only, biallelic, autosomal variants mapped to the GRCh37.p13 (hg19) reference genome.

All analyses were implemented in a reproducible workflow using Snakemake, with parameterization controlled via configuration files to ensure transparency and reproducibility.


## Variant-level quality control

A common variant-level quality control (QC) procedure was applied to all datasets prior to downstream analyses. Only autosomal, SNP-only, biallelic variants were retained, while monomorphic, multiallelic, and palindromic SNPs (A/T and C/G) were removed.

Missingness filters were applied at both variant and sample levels (geno ≤ 0.05; mind ≤ 0.05).

Hardy–Weinberg equilibrium (HWE) filtering was performed separately for each chromosome using a Bonferroni-corrected significance threshold (α = 0.05).


## Relatedness and sample filtering

Sample-level quality control was performed to identify duplicates and related individuals.

Two complementary approaches were used:

- KING, applied without LD pruning and without MAF filtering  
- PLINK IBS estimation (--genome), applied on LD-pruned variants (window = 50 SNPs, step = 5 SNPs, r² = 0.2) with MAF ≥ 0.05  

Based on these analyses, a filtered dataset of unrelated individuals (MASTER_UNRELATED) was defined and used for downstream analyses.


## Population structure analyses

Population structure was investigated using the MASTER_UNRELATED dataset.

### Principal Component Analysis (PCA)

PCA was performed using smartpca (EIGENSOFT) after LD pruning (50 SNP window, 5 SNP step, r² = 0.2) and MAF filtering (MAF ≥ 0.05).

Principal components were used to identify population clustering patterns and assess genetic affinity between the 5Terre population and European reference populations from the 1000 Genomes Project.


### Runs of Homozygosity (ROH)

ROH were computed using PLINK on SNP-only data with MAF ≥ 0.05 and LD pruning optimized for ROH detection (50 SNP window, 5 SNP step, r² = 0.5).

ROH segments were used to infer individual autozygosity burden and were classified based on their length distribution to distinguish between recent and ancient inbreeding signals.


## Linkage disequilibrium decay

LD decay was estimated using combined datasets (5Terre and 1000 Genomes EUR) to determine the genomic scale of correlation between variants.

The decay of LD was evaluated as a function of physical distance, and the characteristic distance at which LD reached one-third of its initial value was estimated. This value was used to define the window size for downstream selection analyses.


## Phasing

Haplotype phasing was performed on the 5Terre dataset using SHAPEIT2. Reference haplotypes from selected European populations (CEU, TSI, IBS) from the 1000 Genomes Project were used as a reference panel.

Only SNPs with MAF ≥ 0.05 were retained for phasing.


## Detection of selection signals (LASSI)

Selection scans were performed using LASSI, a likelihood-based method for detecting selective sweeps from haplotype data.

Input data consisted of phased haplotypes from the 5Terre population. The likelihood window size was defined based on LD decay estimates, and a sliding window approach was applied across the genome.

The method identifies candidate regions under selection and assigns a likelihood-based score to each genomic window. It distinguishes between hard and soft selective sweeps based on the inferred m parameter.


## Post-processing and annotation

Candidate regions identified by LASSI were further processed and annotated through aggregation of overlapping windows, ranking of signals based on likelihood scores, and classification of regions as genic or intergenic.

Selective sweeps were classified as:

- Hard sweeps: m = 1  
- Soft sweeps: m > 1  

Genes overlapping candidate regions were annotated for downstream interpretation.


## Functional enrichment analysis

Genes overlapping candidate regions were analyzed using g:Profiler and STRING to perform pathway enrichment and protein–protein interaction network analysis, providing biological interpretation of the detected signals.


## Optional demographic analysis

An optional analysis of recent demographic history can be performed using IBDNe.

This requires phased genotype data and detection of IBD segments between individuals. Longer IBD segments reflect recent shared ancestry, whereas shorter segments capture more ancient demographic events.

IBDNe estimates effective population size (Ne) over time based on the distribution of IBD segment lengths.





## Reproducibility

All analyses are implemented within a Snakemake workflow to ensure full reproducibility.  
Software environments are managed using Conda, and all parameters are centrally defined in a configuration file.
