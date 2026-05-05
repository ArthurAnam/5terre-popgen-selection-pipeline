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
...
