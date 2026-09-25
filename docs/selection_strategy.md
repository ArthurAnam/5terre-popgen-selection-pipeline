# Selection strategy

## Scope

The selection branch is focused on the Cinque Terre population, with CEU, TSI and IBS used as production comparison populations. FIN and GBR are retained only for European LD context.

The primary statistic is **saltiLASSI Lambda**. Original **LASSI T** is secondary and is used for continuity/benchmarking.

Candidate regions are interpreted as empirical genomic outliers compatible with selective processes, not as definitive adaptive loci or causal variants.

---

## Production marker panels

Production CT/CEU/TSI/IBS selection panels are:

* autosomal
* biallelic SNPs
* phased
* unpruned
* MAF>=0.05 within the population
* A/T and C/G strand-ambiguous SNPs excluded

CT already satisfies the strand-ambiguity rule through core QC. CEU/TSI/IBS are being recalibrated under the same rule so that the marker policy used for LD-window calibration matches the marker policy used for the scan.

A separate CEU Harris & DeGiorgio (2020) benchmark is allowed to use a paper-matched ascertainment distinct from the production panel. It is a methodological benchmark, not the production CEU comparator.

---

## Window calibration

The primary physical LD scale is the first 1-kb bin whose mean hard-call r2 is below one third of the mean r2 for pairs separated by 0.5-1.5 kb.

The five-consecutive-bin crossing is a stability diagnostic only.

For CT:

* pre-phasing MAF>=0.05 panel: 5,007,326 SNPs
* primary LD scale: 55.5 kb
* pre-phasing median at that scale: 99 SNPs
* exact phased scan panel: 4,914,283 SNPs
* exact phased median at 55.5 kb: 97 SNPs
* current salti step candidate: 49 SNPs

Therefore 99 remains provenance for the original calibration; 97/49 is the current exact-input candidate and must be formally re-frozen before production.

Previously reported CEU 116, TSI 112 and IBS 116 SNP windows predate the production A/T/C/G exclusion and are not production-valid until recalculated.

---

## Candidate threshold and region definition

CT uses an empirical-outlier framework:

* primary threshold: top 1% genome-wide Lambda
* descriptive thresholds: top 0.1% and top 5%
* region definition: merge consecutive windows exceeding the primary threshold

Candidate calling is intentionally independent of gene annotation or biological plausibility. Functional information is added only after population-genetic candidate regions have been defined.

Population comparisons use population-specific rankings and locus concordance rather than a common absolute Lambda threshold.

---

## Signal robustness

For each candidate region, inspect:

* maximum and mean Lambda
* local score profile
* genomic span
* number of contributing windows
* local recombination context
* mappability/accessibility
* concordance with secondary LASSI T where informative

The number of overlapping windows is not treated as an independent measure of statistical support because neighboring sliding windows are partially non-independent.

Low-recombination regions receive additional scrutiny because extended LD can broaden composite-likelihood signals.

---

## CRG100 mappability

The production resource is the GRCh37/hg19 UCSC ENCODE CRG100 100-mer alignability BigWig (`wgEncodeCrgMapabilityAlign100mer.bigWig`), whose checksum/resource gate has passed.

The production threshold is **mean CRG100 >=0.9**.

Historical exploratory analyses used `CRG_MIN=0.80` for HIGHCONF and separately defined `CRG_STRONG=0.90`; those values are provenance only.

Before the production filter is encoded, the exact operational unit used for the mean CRG100 calculation (for example the precise scan-window/core interval) must be verified against the saltiLASSI empirical procedure. The threshold is frozen; the interval/application point is not yet frozen.

---

## Functional interpretation

Functional annotation and enrichment are exploratory biological contextualization only. They do not enter candidate calling and are not independent validation of selection.
