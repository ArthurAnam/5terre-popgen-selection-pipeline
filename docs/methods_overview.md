# Methods overview

## Overview

This document summarizes the methodological design of the Cinque Terre population genomics and selection pipeline.

The workflow is designed to support two main analytical components:

* population structure and genomic diversity analyses
* haplotype-based selection scan using LASSI

An optional demographic component based on IBD segment detection and IBDNe may also be performed, depending on data suitability and final analytical priorities.

The selection branch is framed conservatively. Its aim is not to claim definitive adaptive loci or causal variants, but to identify candidate genomic regions showing deviations from neutral expectations and therefore suitable for downstream biological interpretation.

---

## Target population

The target population of the study is the Cinque Terre population.

The expected target dataset consists of:

* whole-genome sequencing data
* 46 WGS samples
* mean depth approximately 30x
* variants already called upstream
* variants already VQSR-filtered upstream

Variant calling and VQSR are considered part of the upstream generation of the Cinque Terre dataset and are not repeated within this workflow.

For this reason, this pipeline does not repeat variant calling, VQSR, depth-based filtering, or genotype-quality-based filtering by default, unless future inspection of the input data suggests that additional filtering is required.

---

## Reference and comparison datasets

The main comparative framework is expected to involve Italian WGS reference cohorts, ideally representing:

* North Italy
* South Italy

These Italian WGS cohorts are currently under evaluation and may be based on previously published Italian genome data, including Sazzini et al. (2020), if available and suitable for integration.

The 1000 Genomes Project EUR panel is also used as a reference resource, mainly for:

* European genetic context
* phasing support
* comparison with selected European populations

The full 1000 Genomes EUR panel includes:

* CEU
* FIN
* GBR
* IBS
* TSI

For selection-related comparison or reference interpretation, the final subset of external populations remains under evaluation.

Candidate 1000 Genomes EUR populations currently include:

* TSI
* FIN
* IBS

The selection scan itself is performed on the Cinque Terre dataset.

External populations are not treated as primary scan targets unless explicitly stated in a future analysis version.

---

## Genome build and coordinate system

All analyses are based on:

* GRCh37.p13 / hg19 coordinates
* VCF 1-based coordinates
* numeric autosomal chromosome naming
* chromosomes 1-22 only

The preferred variant identifier format is:

```text
CHR:POS:REF:ALT
```

This representation is used to reduce ambiguity during dataset harmonization and cross-dataset comparisons.

---

## Dataset harmonization

Before downstream analyses, datasets must be harmonized to ensure that variants are represented consistently across target and reference cohorts.

The harmonization strategy includes:

* using the same genome build
* restricting analyses to autosomal chromosomes
* retaining SNPs only
* excluding multiallelic variants
* excluding monomorphic SNPs
* excluding strand-ambiguous SNPs
* using a consistent variant ID format
* checking allele representation across datasets

This harmonization step is especially important when combining the Cinque Terre dataset with external Italian WGS cohorts and 1000 Genomes EUR reference populations.

---

## Variant and sample filtering strategy

The pipeline applies population-genomics-oriented filtering to the input variant datasets.

The filtering strategy is expressed in terms of removed data:

* remove non-autosomal variants
* remove non-SNP variants
* remove multiallelic variants
* remove monomorphic SNPs
* remove strand-ambiguous SNPs (A/T and C/G)
* remove variants with more than 5% missing genotypes
* remove individuals with more than 5% missing genotypes
* remove variants showing significant deviation from Hardy-Weinberg equilibrium after correction for multiple testing

The current missingness thresholds are:

```text
Variant missingness: 5%
Individual missingness: 5%
```

These thresholds are intended to remove variants and samples with excessive missingness while avoiding overly aggressive filtering in a small WGS dataset.

---

## Hardy-Weinberg equilibrium filtering

Hardy-Weinberg equilibrium filtering is applied chromosome-wise using Bonferroni correction.

For each chromosome, the HWE significance threshold is defined as:

```text
0.05 / number of tested variants on that chromosome
```

This means that the HWE threshold is not a fixed genome-wide value.

Instead, the threshold is computed separately for each chromosome according to the number of variants tested on that chromosome.

This strategy avoids using a single arbitrary HWE threshold across the whole genome and makes the correction explicitly dependent on the number of tests performed within each chromosome.

---

## Relatedness and IBS analyses

Pairwise relatedness and IBS analyses are used mainly for population structure and relatedness exploration.

They also provide a secondary QC role by checking for unexpected duplicates or close relatives.

The current strategy includes:

* KING
* PLINK IBS / `--genome`

The datasets are expected to be composed of unrelated individuals, but this assumption will be checked during the analysis.

KING is used for relatedness estimation without LD pruning and without MAF filtering.

PLINK `--genome` may be used as an additional pairwise IBS check on LD-pruned common variants.

In this workflow, relatedness analyses are not treated as the main QC backbone. Their primary role is exploratory, with a secondary role in confirming that no unexpected close relationships or duplicate samples are present.

---

## Population structure analyses

Population structure analyses are used to place the Cinque Terre population within a broader Italian and European genetic context.

The main comparison is expected to involve:

* Cinque Terre
* North Italy WGS reference cohort
* South Italy WGS reference cohort

The 1000 Genomes EUR panel may provide additional European reference context.

Main analyses include:

* PCA
* runs of homozygosity
* fROH
* total ROH burden
* ROH length distribution
* IBS / pairwise relatedness exploration

These analyses are intended to describe genetic structure, autozygosity, and relationships between Cinque Terre and comparison populations.

The following analyses are not part of the current workflow unless explicitly added in a future version:

* ADMIXTURE
* FST
* ancestry proportion inference
* formal population grouping inference

---

## PCA

Principal component analysis is used to describe genetic structure and to place Cinque Terre individuals relative to Italian and European reference populations.

The current PCA strategy includes:

* SNP-level QC
* MAF filtering
* LD pruning
* PCA using smartpca / EIGENSOFT

Current working parameters:

```text
MAF threshold: 0.05
LD pruning: enabled
Pruning window: 50 SNPs
Pruning step: 5 SNPs
r2 threshold: 0.2
Number of components: 20
```

The PCA is intended as a population structure analysis, not as a formal test of ancestry proportions.

---

## Runs of homozygosity

Runs of homozygosity analyses are used to describe autozygosity and genomic patterns of homozygosity in the Cinque Terre population.

Main outputs include:

* ROH segments
* fROH
* total ROH burden
* ROH length distribution

ROH analyses are useful for comparing the distribution and burden of homozygous segments between Cinque Terre and external reference populations.

The current working strategy includes:

* MAF filtering
* LD pruning
* PLINK `--homozyg`

Current working parameters:

```text
MAF threshold: 0.05
LD pruning: enabled
Pruning window: 50 SNPs
Pruning step: 5 SNPs
r2 threshold: 0.5
```

Final PLINK `--homozyg` parameters remain to be refined based on SNP density, WGS data characteristics, and comparison dataset harmonization.

---

## Optional demographic analysis

Recent demographic history may be explored using IBD-based approaches.

This branch is optional and depends on data suitability.

Potential workflow:

* phased haplotypes
* hap-IBD or similar IBD segment detection
* IBDNe for recent effective population size inference

IBDNe estimates are interpreted in generations and may be converted into calendar years using a generation time of 29 years.

This conversion affects only the temporal interpretation of the results, not the genetic inference itself.

The optional demographic branch is not required for the main population structure and selection analyses.

---

## Selection scan

The selection branch is focused on the Cinque Terre population.

The main selection method is:

* LASSI

The selection scan is based on phased haplotype data and is used to identify candidate genomic regions showing deviations from neutral expectations.

Main steps include:

* MAF filtering
* LD decay estimation
* phasing with SHAPEIT2
* LASSI genome-wide scan using a sliding likelihood window
* candidate region definition
* candidate region ranking
* functional genomic annotation

External populations may be used for comparison, reference interpretation, or phasing support, but the primary LASSI scan target is the Cinque Terre dataset.

Candidate external populations for selection-related comparison currently include:

* TSI
* FIN
* IBS

The final comparison set remains under evaluation.

---

## LD decay and LASSI window definition

LD decay is used to derive suitable LASSI window-size parameters.

The current working strategy estimates LD decay and uses the decay pattern to define the SNP-based window length for LASSI.

Current working parameters:

```text
Baseline distance: 1 kb
Decay fraction: one third of baseline LD
Shift fraction: approximately 10% of the LASSI window
```

The derived parameters are:

* LASSI window size in SNPs
* LASSI shift size in SNPs

These values are dataset-dependent and therefore remain to be finalized after LD decay estimation.

LD decay and phasing are both preparatory steps for the selection branch.

Phasing is not considered a downstream consequence of LD decay.

---

## Phasing

The Cinque Terre dataset is phased before LASSI analysis.

The current phasing method is:

* SHAPEIT2

The 1000 Genomes Project EUR panel is expected to provide primary phasing/reference support.

Additional Italian WGS reference cohorts may also be evaluated as potential support, depending on data availability, compatibility, and final analytical design.

Phasing parameters remain to be refined, including:

* effective population size
* burn-in iterations
* pruning iterations
* main iterations
* genetic map

---

## LASSI

LASSI is used as a likelihood-based method to scan phased haplotype data for candidate genomic regions showing deviations from neutral expectations.

The LASSI scan is interpreted as a genome-wide sliding-window likelihood scan.

The primary output is not a definitive list of causal adaptive loci, but a set of candidate regions requiring downstream interpretation.

Current working outputs include:

* LASSI score profiles
* high-scoring windows
* candidate genomic regions
* region-level summaries
* functional genomic annotations

The LASSI window size and shift size are derived from LD decay analysis and remain to be finalized.

---

## Candidate region definition

Candidate regions may consist of:

* a single high-scoring window
* multiple nearby high-scoring windows merged into a candidate region

Candidate region merging must consider:

* LASSI score magnitude
* local signal coherence
* maximum distance in kb allowed between neighboring windows
* genomic span
* functional genomic context

The genomic span of LASSI windows and candidate regions is not directly provided by LASSI.

It must be computed post hoc from genomic coordinates.

Candidate regions are interpreted as genomic regions compatible with possible selective processes, not as definitive proof of adaptive loci or causal variants.

---

## Candidate region ranking

Candidate regions should be ranked using region-level summaries rather than individual windows alone.

Possible ranking criteria include:

* maximum LASSI score
* mean LASSI score
* local signal coherence
* genomic span
* overlap with functional genomic elements
* biological interpretability

The final ranking strategy remains to be refined once real LASSI outputs are available.

---

## Functional genomic annotation

Candidate regions may be annotated according to their functional genomic context.

The main functional classes are:

* coding regions
* regulatory regions
* intergenic regions

This annotation is used to contextualize candidate regions and support biological interpretation.

Functional genomic context should not be interpreted as independent evidence of selection.

---

## Exploratory biological contextualization

Functional interpretation is treated as exploratory biological contextualization.

ORA and/or GSEA-like approaches may be used to provide biological context for genes overlapping or near candidate regions.

These analyses are not considered independent validation of selection.

They are used to support interpretation and hypothesis generation.

The main result of the selection branch remains the definition and prioritization of candidate genomic regions deviating from neutral expectations.

---

## Reproducibility principles

The workflow is designed to be modular and reproducible.

The main reproducibility principles are:

* central configuration in `config/config.yaml`
* explicit documentation of parameters and assumptions
* separation between raw data, processed data, and results
* use of Snakemake for workflow execution
* use of Conda environments for software dependencies
* preservation of logs and benchmark outputs
* avoidance of manual, undocumented intermediate file editing

Raw data should not be modified directly.

Generated outputs should be reproducible from input data, configuration files, environment files, and workflow rules.

---

## Output organization

The expected output structure is:

```text
results/
logs/
benchmarks/
```

The `results/` directory stores generated outputs.

The `logs/` directory stores command-level logs.

The `benchmarks/` directory stores runtime and resource usage information when available.

Raw input data should remain outside generated output directories and should not be overwritten by the workflow.

---

## Notes on interpretation

The pipeline is designed to support careful population-genetic interpretation.

Population structure analyses are descriptive and comparative.

Selection scans identify candidate genomic regions deviating from neutral expectations.

Functional interpretation provides exploratory biological context.

No single downstream analysis is treated as definitive proof of adaptation or causality.

Interpretation should therefore integrate:

* statistical signal strength
* genomic context
* comparison with reference populations
* local signal coherence
* biological plausibility
* technical reliability of the region
