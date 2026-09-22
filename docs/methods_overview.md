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

## Relatedness analysis

Pairwise relatedness is evaluated within the Cinque Terre cohort using KING.
This analysis has a secondary QC role by checking for unexpected duplicates or
close relatives before downstream population analyses.

The current strategy is:

* KING-Robust pairwise kinship using `--kinship`
* no KING-specific MAF filter
* no LD pruning
* inspection of the kinship coefficient together with IBS0
* no automatic sample exclusion

KING is run on the stable target-QC dataset before final HWE filtering. This is
a project-specific ordering choice rather than an explicit KING requirement:
if relatedness leads to a sample exclusion, cohort-dependent site statistics
and HWE can then be recalculated on the final retained cohort.

Duplicate/MZ, first-degree and second-degree relationships are the primary
screening targets. The standard 0.0442 third-degree boundary is retained only
as an exploratory flag because the KING documentation describes `--kinship`
as most reliable for closer relationships.

`--related` and `--ibdseg` are reserved as follow-up analyses if the primary
kinship screen identifies a close or ambiguous pair requiring IBD-based
refinement.

The supporting paper and software documentation are recorded in
`docs/methodological_references.md`.

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

For the harmonized Cinque Terre + reference PCA dataset, the primary marker
filter is fixed at MAF >= 0.05 followed by LD pruning with 50-SNP windows,
5-SNP steps and an r^2 threshold of 0.2. This matches the PCA preprocessing
used by Sazzini et al. (2020) in a closely related Italian population-genomics
study. The number and genomic distribution of retained markers will be recorded
after pruning. The number of PCs and the final joint-PCA versus projection
strategy remain to be finalized after inspecting the harmonized dataset.

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

The production method is PLINK `--homozyg` on the joint CT + 1000G EUR
MAF >= 0.05 marker panel without LD pruning. The frozen call uses a 50-SNP
sliding window, a minimum of 50 SNP per called ROH, a minimum physical length
of 1.5 Mb, density <=50 kb/SNP, gap <=500 kb, <=5 missing calls per window,
<=1 heterozygous call per window, and window hit threshold 0.05.

Primary fROH is the summed length of autosomal ROH >=1.5 Mb divided by
2.77e9 bp. Light-VIF pruning, gap/heterozygote settings, and marker-frequency
scope were evaluated in prespecified sensitivity analyses before freezing the
production definition.

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

The production LASSI comparison set is:

* Cinque Terre (focal discovery population)
* CEU (empirical LASSI benchmark)
* TSI (Southern-European comparator)
* IBS (Southern-European comparator)

FIN and GBR are retained for European LD context but are not part of the primary LASSI scan.

---

## LD decay and LASSI window definition

LD decay is used to calibrate the physical scale of the SNP-delimited LASSI
window. The procedure follows the empirical human protocol of Harris and
DeGiorgio (2020): LD is measured as pairwise r^2 and the relevant physical
interval is the first distance at which mean LD falls below one third of the
value observed for SNP pairs separated by approximately 1 kb.

Because the final CT WGS dataset is extremely dense, an all-pairs LD
calculation would be unnecessarily large. The audit therefore selects 100
anchor SNPs per autosome at deterministic, approximately even physical
positions (excluding 500-kb chromosome edges) and computes unphased hard-call
r^2 between each anchor and all SNPs within 500 kb. LD is summarized in
1-kb distance bins. The primary calibration uses MAF >=0.05 and a prespecified
MAF >=0.01 sensitivity is run in parallel.

LD-decay calibration is complete. Population-specific physical decay scales
were converted to empirical SNP counts in the corresponding MAF>=0.05 panels.
The frozen LASSI winsize/winstep values are CT 99/10, CEU 116/12, TSI 112/11,
and IBS 116/12 SNPs. FIN and GBR are retained as LD-context-only populations.

LD decay and phasing are both preparatory steps for the selection branch;
phasing is not a downstream consequence of LD decay.

---

## Phasing

The Cinque Terre dataset is phased before LASSI analysis.

The production phasing software is SHAPEIT5 `phase_common`, replacing the
obsolete SHAPEIT2 plan. SHAPEIT5 is the maintained successor and is designed
for common-variant phasing as the first stage of WGS phasing.

Before production phasing, the workflow quantifies exact CHR:POS:REF:ALT
overlap between the CT MAF>=0.05 panel and two possible phased reference scopes:
the full 1000 Genomes Phase 3 panel and the 503-sample EUR subset. This
preflight is required because SHAPEIT5 reference-assisted phasing considers
only target variants present in the reference panel.

The provisional production design is two-stage: first phase exact
target-reference shared variants with reference support; then use the resulting
CT haplotypes as a scaffold to phase the complete CT MAF>=0.05 panel. This
preserves target common variants that are absent from the external reference
while still exploiting external haplotype information.

GRCh37 SHAPEIT5 genetic maps will be used. Reference scope and whether
whole-chromosome or chunked phasing is preferable remain to be finalized after
the overlap preflight.

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

The LASSI window size and shift size are derived from population-specific LD-decay calibration and are frozen at CT 99/10, CEU 116/12, TSI 112/11 and IBS 116/12 SNPs.

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
