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

The target population is the Cinque Terre WGS cohort. The delivered dataset contains 50 individuals at approximately 30x mean depth. Variant calling, joint genotyping, restriction to PASS biallelic SNPs and VQSR were completed upstream.

The current workflow reconstructs sample-level QC from all 50 delivered individuals and documents the exclusion of four reviewed outliers. The retained downstream cohort therefore contains **46 individuals**.

No additional arbitrary site/genotype QUAL, DP, GQ or allele-balance hard thresholds are imposed after those reviewed exclusions. These annotations were used diagnostically during sample review, while downstream site QC is based on the population-genomic criteria documented below.

---

## Reference and comparison datasets

The implemented European reference panel is the 1000 Genomes Project Phase 3 EUR subset: CEU (99), FIN (99), GBR (91), IBS (107) and TSI (107), for 503 individuals in total.

The panel serves three distinct roles:

* European context for population-structure and ROH analyses;
* EUR reference support for CT SHAPEIT2 phasing;
* population-specific selection comparison.

The production selection comparison set is **CEU, TSI and IBS**. FIN and GBR are retained as LD-context populations and are not primary selection-scan comparators.

Production CEU/TSI/IBS selection panels are processed population-wise with the same common-variant/strand-ambiguity policy used for CT. A separate CEU methodological benchmark will be kept distinct from the production CEU scan and will be configured to approximate the empirical Harris & DeGiorgio (2020) LASSI-T analysis as closely as the available data allow.

Italian WGS cohorts remain a possible future extension for broader population-structure comparisons; they are not required for the current production selection branch.

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

Cross-dataset analyses use GRCh37/hg19 coordinates, autosomes 1-22 and the variant identifier `CHR:POS:REF:ALT`.

For CT + 1000G EUR harmonization, the EUR subset is restricted to biallelic polymorphic SNPs, sites with >5% missingness are excluded, and strand-ambiguous A/T and C/G SNPs are removed before exact allele matching. Only exact `CHR:POS:REF:ALT` matches are retained; no strand-flip or allele-rescue procedure is used.

The delivered CT input is already a PASS biallelic-SNP dataset from the upstream calling workflow. Therefore “SNP/biallelic filtering” is treated as an input validation for CT rather than a newly imposed downstream filter.

---

## Variant and sample filtering strategy

After the four reviewed sample exclusions, CT core QC is:

* autosomes 1-22;
* removal of A/T and C/G strand-ambiguous SNPs;
* removal of sites that are monomorphic in the retained cohort;
* removal of sites with >5% missing genotypes;
* evaluation of individual missingness at >5%;
* repetition of cohort-dependent site QC after any individual removal until stable;
* chromosome-wise HWE filtering after the stable cohort is obtained.

No general MAF threshold is part of core WGS QC. MAF thresholds are analysis-specific downstream choices.

Because the delivered CT VCF is already composed of PASS biallelic SNPs, non-SNP/multiallelic checks validate the delivered input rather than constituting additional filtering introduced by this workflow.

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

The implemented reference framework for the current analyses is CT + 1000 Genomes EUR (549 individuals after CT QC).

The main completed analyses are:

* joint PCA;
* KING relatedness within CT;
* runs of homozygosity and fROH;
* ROH burden and length-distribution comparisons;
* exploratory pairwise population differentiation used only as context.

Italian WGS cohorts can be added later if suitable data become available, but they are not required to interpret the currently completed CT + EUR analyses.

ADMIXTURE and ancestry-proportion inference are not part of the current workflow.

---

## PCA

PCA is performed jointly on CT and the five 1000 Genomes EUR populations using smartpca/EIGENSOFT.

The joint MAF>=0.05 panel contains 4,878,327 SNPs before LD pruning. The primary pruning procedure uses 50-SNP windows, 5-SNP steps and r2=0.2, retaining 295,375 markers. A sensitivity analysis additionally masks canonical high-LD regions and retains 287,815 pruned markers.

The first three PCs are essentially unchanged by the high-LD-region sensitivity, supporting their use for the main population-structure interpretation. PCA remains descriptive; it is not interpreted as a formal ancestry-proportion estimate.

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

The primary selection method is **saltiLASSI Lambda**, implemented with the maintained `lassip` v1.2.1 codebase. Original **LASSI T** is retained as a secondary benchmark/continuity statistic.

The focal discovery population is CT. CEU, TSI and IBS are the production comparison populations; FIN and GBR remain European LD-context populations only.

Production marker panels are population-specific, phased, unpruned and use MAF>=0.05. Strand-ambiguous A/T and C/G SNPs are excluded from the production CT/CEU/TSI/IBS selection panels. CEU additionally has a separate paper-oriented LASSI-T benchmark branch, which must not be confused with the production matched CEU comparator.

Candidate calls are based on within-population score distributions rather than absolute score equality across populations.

---

## LD decay and LASSI window definition

LD decay calibrates the physical scale that is then translated into a SNP-delimited local HFS window. For each population, unphased hard-call r2 is summarized from deterministic physical anchors. The primary distance is the first 1-kb bin after the 0.5-1.5-kb baseline whose mean r2 falls below one third of the baseline mean. A five-consecutive-bin crossing is retained only as a stability diagnostic.

For CT, the MAF>=0.05 pre-phasing panel gave a primary scale of 55.5 kb and a median of 99 SNPs. SHAPEIT2 reference alignment subsequently reduced the exact phased scan panel from 5,007,326 to 4,914,283 SNPs. Recounting the same 2,200 anchors on that exact phased panel gives a median of **97 SNPs** at 55.5 kb (49-SNP approximately half-window step candidate); the stable 57.5-kb diagnostic gives a median of 101 SNPs.

Accordingly, 99 is retained as the pre-phasing calibration result, while **97/49 is the current exact-phased-panel candidate pending formal production re-freeze**.

The previously reported CEU 116, TSI 112 and IBS 116 SNP windows were calibrated before the production-wide decision to exclude A/T and C/G SNPs from the selection comparators. Those values are historical/pre-filter calibration results and must be recalculated on the final non-palindromic CEU/TSI/IBS panels before production scanning.

LD decay and phasing are independent preparatory steps; phasing is not a consequence of the LD calculation.

---

## Phasing

The Cinque Terre dataset is phased before LASSI using SHAPEIT2 (v2.r904) with
the phased 1000 Genomes Phase 3 EUR reference panel.

The CT cohort contains 46 individuals. The existing local reference files contain all 2503 Phase 3 individuals; the accompanying SAMPLE file labels each individual by population and super-population, enabling SHAPEIT2 to restrict the reference internally to EUR without generating a second physical HAP/LEGEND panel. SHAPEIT2 documentation specifically
states that reference-assisted phasing is particularly useful for studies with
fewer than approximately 100 individuals. The European/CEU effective
population-size value recommended by SHAPEIT2, Ne=11,418, is therefore used.

The target phasing panel contains biallelic CT SNPs with MAF>=0.05, matching
the common-variant panel used for CT LASSI. A chromosome-wise overlap preflight is followed by a formal two-stage SHAPEIT2 `-check` audit against the EUR reference before phasing. The completed audit retained 4,914,283 of 5,007,326 CT MAF>=0.05 SNPs (98.14%) after excluding 90,392 study variants absent from the EUR reference and 2,651 allele-misaligned variants; all 22 post-exclusion checks passed.
Study SNPs absent from the reference or showing incompatible alleles are
excluded from the reference-assisted run according to the SHAPEIT2 alignment
output.

For sequence-derived genotypes, the SHAPEIT2 documentation recommends a
0.5-Mb conditioning window rather than the 2-Mb default used for typical GWAS
data, based on the developers' sequencing experiments. This choice is also
supported by Delaneau et al. (2013, *American Journal of Human Genetics*),
who used the standard SHAPEIT2 model with W=0.5 Mb when phasing high-coverage
sequence genotypes together with a European 1000 Genomes reference, and by
Sharp et al. (2016, *Bioinformatics*), who state that SHAPEIT2 had previously
shown good performance at 0.5 Mb for unphased genotypes derived from sequencing
and therefore used that window in their reference-based sequencing experiments.
The production candidate therefore uses `--window 0.5`,
`--effective-size 11418`, 400 conditioning states, and the documented
default MCMC schedule (7 burn-in, 8 pruning and 20 main iterations). The
400-state setting was frozen after a chromosome-20 benchmark confirmed
practical runtime and sub-1-GB peak resident memory with one thread. The `--no-mcmc` shortcut is not used because SHAPEIT2 recommends
it only for much smaller study samples, typically fewer than 10 individuals.

A fixed random seed and one SHAPEIT2 thread per chromosome are used for exact
reproducibility; chromosome jobs may be parallelized by Snakemake. GRCh37
genetic maps are supplied explicitly.

CEU, TSI and IBS are already phased in the 1000 Genomes source and are not
re-phased for LASSI.

---

## saltiLASSI and original LASSI

The primary scan is saltiLASSI Lambda. The maintained `lassip` v1.2.1 implementation is pinned to commit `a6a9d18c2323330fbf74d5a490f9e9c4ebe41d7c`.

Core model settings are phased input, K=10, Model D with sweeping-class mass proportional to exp(-i) (`--lassi-choice 4` in lassip v1.2.1), physical-distance mode and `--max-extend-bp 100000`.

Original LASSI T is run secondarily. CEU has two explicitly different roles: a production matched comparator processed under the same marker policy as CT/TSI/IBS, and a separate Harris & DeGiorgio (2020)-oriented T-statistic benchmark. Results from those branches must not be pooled.

Absolute Lambda/T values are not assumed to be directly comparable across populations with different haplotype-frequency spectra or ascertainment; comparative interpretation emphasizes within-population ranks, candidate loci and concordance.

---

## Candidate region definition

CT uses a pre-specified empirical-outlier framework because a sufficiently calibrated CT demographic model is not available for a defensible demography-matched whole-genome neutral threshold.

The primary CT threshold is the top 1% of genome-wide Lambda. Top 0.1% and top 5% are descriptive summaries. Consecutive windows above the primary threshold are merged into candidate/outlier regions.

These regions are called **candidate/outlier regions**, not simulation-based significant regions. Region calling is based on the score threshold and window adjacency; functional annotation is applied only after candidate definition and is not allowed to determine whether a region is called.

---

## Candidate region ranking

Candidate regions are summarized using population-genetic quantities such as maximum/mean score, genomic span, number of contributing windows and local score profile. Low-recombination or technically problematic regions receive additional scrutiny.

Functional annotations may be attached for interpretation but are not used as independent evidence that a region is under selection.

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


The production workflow therefore requires a clean second SHAPEIT2 `-check`
after applying the diagnostic exclusion list. A non-zero first check is not
treated as a pipeline failure when SHAPEIT2 has generated both the detailed
`.snp.strand` report and the corresponding `.snp.strand.exclude` file;
any other non-zero termination remains fatal. The phasing rule depends on the
per-chromosome post-exclusion PASS sentinel, preventing phasing from starting
unless the reference alignment has been revalidated.
