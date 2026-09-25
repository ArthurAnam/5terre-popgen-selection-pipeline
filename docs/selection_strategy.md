# Selection strategy

## Populations

Production scans:

- CT
- CEU
- TSI
- IBS

LD context only:

- FIN
- GBR

CEU also has a separate Harris & DeGiorgio (2020)-oriented LASSI-T benchmark.

## Marker panels

Production CT/CEU/TSI/IBS panels are:

- autosomal
- biallelic SNPs
- phased
- unpruned
- MAF >=0.05 within population
- A/T and C/G excluded

The LD-window calibration must use the same marker policy as the final scan panel.

## Window calibration

Primary LD scale:

`d_LD,p = min{d: mean_r2_p(d) < (1/3) * mean_r2_p(1 kb)}`

Operationalization:

- 0.5-1.5 kb baseline
- 1-kb bins
- first crossing = primary
- five consecutive crossings = stability diagnostic

CT:

- primary physical scale: 55.5 kb
- pre-phasing calibration: 99 SNPs
- exact phased panel: 4,914,283 SNPs
- exact phased median at 55.5 kb: 97 SNPs
- current saltiLASSI step candidate: 49 SNPs

CEU/TSI/IBS windows are recalibrated after A/T and C/G removal.

## Scan

Primary:

- saltiLASSI Lambda

Secondary:

- original LASSI T

CT candidate threshold:

- top 1% genome-wide Lambda

Descriptive thresholds:

- top 0.1%
- top 5%

Consecutive above-threshold windows are merged into candidate/outlier regions.

## CRG100

Production threshold:

- mean CRG100 >=0.9

Historical exploratory thresholds:

- CRG_MIN=0.80
- CRG_STRONG=0.90

The exact production interval used for the CRG100 mean is still being verified.

## Interpretation

Candidate regions are empirical population-genetic outliers.

Functional annotation is added after candidate definition and is treated as exploratory context, not independent evidence of selection.
