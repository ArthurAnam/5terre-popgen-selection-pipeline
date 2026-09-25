# Pipeline diagram

## Current workflow

```mermaid
flowchart TD
    A[Delivered CT WGS<br/>50 samples<br/>11,490,646 PASS SNPs]
    A --> B[Sample QC review]
    B --> C[46 retained CT samples]

    C --> D[CT core QC<br/>autosomes<br/>remove A/T & C/G<br/>monomorphic + missingness<br/>chr-wise HWE Bonferroni]
    D --> E[9,000,246 CT SNPs]

    R[1000G Phase 3 EUR<br/>503 samples<br/>CEU FIN GBR IBS TSI]

    E --> H[CT + EUR exact harmonization]
    R --> H
    H --> P[PCA / ROH / EUR context]

    E --> M[CT MAF>=0.05<br/>5,007,326 SNPs]
    M --> LDCT[CT LD decay<br/>55.5 kb primary]
    M --> S[SHAPEIT2 reference-assisted phasing]
    R --> S
    S --> SA[Formal check + post-phasing audit<br/>4,914,283 phased SNPs]
    SA --> WD[Exact phased-density audit<br/>97-SNP current candidate at 55.5 kb]

    R --> REF[CEU / TSI / IBS production panels<br/>population MAF>=0.05<br/>remove A/T & C/G]
    REF --> LDR[Recalibrate population-specific LD windows]

    WD --> L[Production saltiLASSI Lambda]
    LDR --> L
    L --> T[Original LASSI T secondary]
    L --> O[CT top-1% empirical outliers<br/>merge consecutive windows]
    O --> Q[CRG100 mean >=0.9<br/>exact interval definition pending verification]
    Q --> F[Functional annotation<br/>exploratory context only]

    R -. CEU paper-matched ascertainment .-> BENCH[Separate H&D2020-oriented<br/>LASSI T benchmark]
```

## Interpretation

The delivered CT dataset is not re-called or re-VQSR-filtered. The workflow begins from the upstream PASS biallelic-SNP VCF and applies population-genomic QC.

The 1000 Genomes EUR panel has different roles in different branches: joint European context, SHAPEIT2 reference support, and population-specific CEU/TSI/IBS selection comparators.

Production selection comparison requires matched marker-policy logic across CT/CEU/TSI/IBS. The separate CEU Harris & DeGiorgio benchmark is intentionally not forced to share the production ascertainment because its purpose is methodological continuity with the published original-LASSI empirical analysis.

Italian WGS cohorts remain a possible future extension and are not a dependency of the current production selection branch.
