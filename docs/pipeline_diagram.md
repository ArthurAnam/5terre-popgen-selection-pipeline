# Pipeline diagram

```mermaid
flowchart TD

A[Raw / Harmonized datasets] --> B[Variant-level QC]
B --> C[MASTER QC DATASET]

C --> D[Relatedness / IBS QC]
D --> E[MASTER_UNRELATED]

E --> F[PCA]
E --> G[ROH]

C --> H[LD decay]
H --> I[LASSI window definition]

C --> J[Selection branch]
J --> K[MAF ≥ 0.05]
K --> L[Phasing]
L --> M[LASSI]
M --> N[Annotation]
N --> O[ORA]

E --> P[Optional: IBD / IBDNe]

```

Description

The pipeline is structured around a central quality-controlled dataset, from which multiple analysis branches originate.

After variant-level quality control, a master dataset is generated and further filtered for relatedness to obtain a set of unrelated individuals. This dataset is used for population structure analyses, including PCA and ROH.

In parallel, the QC-filtered dataset is used for selection analyses. Linkage disequilibrium (LD) decay is estimated to define the genomic window size used in LASSI. The target dataset is then phased and analyzed using a likelihood-based framework to identify candidate regions under selection.

An optional branch based on IBD segment detection and IBDNe can be used to infer recent demographic history.
