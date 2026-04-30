# Pipeline diagram

```text
Raw / harmonized datasets
        ↓
Variant-level QC backbone
        ↓
MASTER QC DATASET
        ↓
├── Relatedness / IBS QC
│       ↓
│   MASTER_UNRELATED
│       ↓
│   ├── PCA-A: 5Terre + 1KG EUR
│   ├── PCA-B: 5Terre + 1KG EUR + Human Origins
│   └── ROH
│
├── LD decay
│       ↓
│   LASSI window parameter
│
└── Selection branch
        ↓
    MAF ≥ 0.05
        ↓
     Phasing
        ↓
      LASSI
        ↓
    Annotation
        ↓
       ORA
