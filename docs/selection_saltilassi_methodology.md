# saltiLASSI pre-production methodology

## Current status

The model-level design is frozen, but the SNP-window component is partially reopened while the final production marker panels are audited.

Primary method: **saltiLASSI Lambda**  
Secondary method: **original LASSI T**  
Implementation: `lassip` v1.2.1, commit `a6a9d18c2323330fbf74d5a490f9e9c4ebe41d7c`

No production genome-wide saltiLASSI scan has started.

---

## Core model settings

* phased input
* K=10
* Model D, sweeping-class mass proportional to exp(-i)
* `--lassi-choice 4` in lassip v1.2.1
* physical distance (`bp`)
* `--max-extend-bp 100000`
* population-specific null HFS
* null pK floor for K=10: 0.001

The original Harris & DeGiorgio Python option numbering must not be transferred directly to lassip: the paper's Model D corresponds to `--lassi-choice 4` in the maintained v1.2.1 implementation.

---

## Production marker ascertainment

Production CT/CEU/TSI/IBS panels use MAF>=0.05 within population and exclude A/T and C/G strand-ambiguous SNPs. Panels are not LD-pruned.

This differs from the 2022 empirical CEU/YRI saltiLASSI analysis, which used sample-polymorphic biallelic SNPs. Therefore published CEU/YRI absolute Lambda values are not treated as directly comparable to the production scans.

CEU also has a separate Harris & DeGiorgio (2020)-oriented original-LASSI-T benchmark. That benchmark may use a paper-matched ascertainment distinct from the production CEU panel and must be labelled separately.

---

## LD-to-SNP window calibration

For population p:

`d_LD,p = min{d: mean_r2_p(d) < (1/3) * mean_r2_p(1 kb)}`

The 1-kb reference is operationalized as pairs separated by 0.5-1.5 kb. LD is summarized in 1-kb bins. The first below-threshold bin is primary; a five-consecutive-bin crossing is a robustness diagnostic only.

### CT

The original MAF>=0.05 pre-phasing panel contained 5,007,326 SNPs. Its 55.5-kb primary LD scale corresponded to a genome-wide median of 99 SNPs.

After the formal SHAPEIT2 alignment check, the exact phased scan panel contains 4,914,283 SNPs. The lassip-input audit passes on all 22 autosomes, 46 samples, zero missing genotypes, zero unphased called genotypes, zero non-binary genotypes, zero duplicate sites and zero ordering failures.

Recounting the original 2,200 physical anchors on the exact phased panel gives:

* 55.5-kb primary width: median 97 SNPs, mean 101.55, IQR 66-134
* current approximately-50% salti step candidate: 49 SNPs
* 57.5-kb stability width: median 101 SNPs

Thus **97/49 is the current exact-phased-panel saltiLASSI candidate**. The older 99/50 value remains provenance until the production value is formally re-frozen.

For secondary LASSI T, the step remains a separate ~10%-of-window design choice and must be updated consistently once the CT production winsize is finalized.

### CEU, TSI and IBS

The previous 116/58, 112/56 and 116/58 salti values were derived before the production decision to exclude A/T and C/G SNPs from the comparator marker panels.

Those values are therefore **pre-filter historical calibration results**, not current production settings. LD decay and SNP-window conversion must be rerun on the final non-palindromic MAF>=0.05 population-specific panels.

---

## Null HFS

If no external null spectrum is supplied, lassip computes a population-specific genome-wide mean K-truncated HFS across supplied windows. Observed and null spectra must use the same marker ascertainment within each population.

Comparisons among CT, CEU, TSI and IBS emphasize population-specific ranks and candidate-locus concordance, not equality of raw Lambda values.

---

## CRG100

The verified production resource is UCSC hg19 `wgEncodeCrgMapabilityAlign100mer.bigWig`; its resource/checksum gate passes.

The production threshold is mean CRG100 >=0.9. Historical exploratory analyses used an operational HIGHCONF cutoff of 0.8 and separately labelled 0.9 as CRG_STRONG.

The remaining methodological item is the exact genomic unit over which the CRG100 mean is calculated and the exact point at which that filter is applied. This must be verified before production coding.

---

## Inference framework

CT does not currently have a sufficiently calibrated demographic model for a defensible demography-matched whole-genome neutral significance threshold.

Therefore:

* primary CT candidate threshold: top 1% genome-wide Lambda
* descriptive thresholds: top 0.1% and top 5%
* merge consecutive above-threshold windows
* terminology: candidate/outlier regions, not simulation-based significant regions

Original LASSI T is secondary. CEU's paper-oriented T benchmark is evaluated separately from the production matched CEU scan.
