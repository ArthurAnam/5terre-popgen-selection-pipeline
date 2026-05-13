# Selection strategy

## Candidate region definition

Selective sweep candidates are identified from genome-wide LASSI likelihood profiles generated on phased haplotype data.

Candidate regions are interpreted by considering both high-scoring windows and local consistency of the surrounding signal, since the true selective target may not necessarily coincide with the single highest-scoring genomic window.

The identification strategy combines:

- likelihood score magnitude
- spatial consistency across neighboring windows
- overlap between adjacent windows
- biological interpretability

Current open methodological decisions include:

- genome-wide vs chromosome-wise thresholding
- empirical percentile thresholds
- definition of candidate region boundaries
- handling of overlapping windows
- maximum distance allowed for candidate region merging
- minimum genomic span required for candidate regions
- ranking strategy for candidate sweeps

The current working strategy is based on empirical ranking of high-scoring windows, with additional consideration of overlap consistency and biological plausibility.

Previous exploratory analyses used:

- top 1% empirical thresholds
- thresholds calibrated from CEU simulations reported in the original LASSI publication
- comparison with populations showing partially similar effective population size decay patterns

These decisions remain subject to refinement during the final analysis phase.

---

## Candidate signal robustness

High-scoring LASSI windows will be evaluated not only by their likelihood score, but also by the spatial stability of the surrounding signal.

Particular attention will be given to the relationship between genomic span and signal intensity. Narrow isolated peaks may reflect stochastic variation or local LD patterns, but they may also represent older or spatially restricted selective events.

For this reason, candidate regions will be assessed using:

- LASSI score magnitude
- genomic span
- maximum distance allowed for candidate region merging
- overlap consistency across adjacent windows
- local signal continuity
- genomic accessibility / mappability
- functional genomic context

The distribution of candidate region widths will be inspected to distinguish broad, spatially coherent signals from isolated narrow peaks.

The number of overlapping windows will not be interpreted as an independent measure of robustness, since it is strongly influenced by SNP density and window construction strategy.

Very narrow candidate regions will not be automatically discarded, since selective signals may differ substantially in genomic extent depending on recombination patterns, sweep age, local genomic architecture, and demographic history.

---

## Genomic accessibility and problematic regions

Candidate regions will be evaluated in relation to genomic accessibility and mappability.

Regions overlapping poorly mappable or otherwise problematic genomic intervals may produce unreliable signals due to mapping uncertainty, phasing errors, local alignment ambiguity, or variant calling artifacts.

Following the strategy adopted in the original LASSI framework, candidate windows overlapping regions of low alignability and mappability may be filtered using CRG100 scores.

Previous exploratory analyses excluded windows overlapping genomic regions with mean CRG100 score < 0.9.

This filtering step is intended to reduce false positive signals arising from technically unreliable genomic regions, rather than to impose biological assumptions about selection.

Problematic regions may therefore be:

- excluded before final candidate prioritization
- retained but flagged as lower-confidence candidates

depending on the final validation strategy.
