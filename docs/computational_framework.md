# Computational framework

## Main software

| Tool | Version / implementation | Main role |
|---|---|---|
| Snakemake | 9.27.0 | workflow orchestration |
| bcftools | pipeline environment | VCF filtering and manipulation |
| PLINK 1.9 / 2 | pipeline environment | QC, LD, PCA preprocessing, ROH |
| KING | 2.3.2 | relatedness |
| EIGENSOFT / smartpca | 7.2.1 | PCA |
| SHAPEIT2 | v2.r904 | CT phasing |
| lassip | v1.2.1, commit `a6a9d18c2323330fbf74d5a490f9e9c4ebe41d7c` | saltiLASSI Lambda and original LASSI T |

## SHAPEIT2

Used for reference-assisted phasing of CT with the 1000 Genomes Phase 3 EUR panel.

Main production settings:

- Ne=11,418
- 0.5-Mb conditioning window
- 400 states
- 7/8/20 burn-in/pruning/main iterations
- fixed seed
- one thread per chromosome

Reference:
Delaneau O, Zagury JF, Marchini J. 2013. *Nature Methods* 10:5-6. https://doi.org/10.1038/nmeth.2307

## lassip

Primary scan:

- saltiLASSI Lambda via `--salti`

Secondary scan:

- original LASSI T via `--lassi`

Production model:

- phased input
- K=10
- Model D
- `--lassi-choice 4`
- physical-distance mode
- `--max-extend-bp 100000`

References:

- Harris AM, DeGiorgio M. 2020. *Molecular Biology and Evolution* 37:3023-3046. https://doi.org/10.1093/molbev/msaa115
- DeGiorgio M, Szpiech ZA. 2022. *PLOS Genetics* 18:e1010134. https://doi.org/10.1371/journal.pgen.1010134

## Other tools

- GNU awk, sed, grep, coreutils, Bash: text and workflow utilities
- tabix: compressed genomic indexing
- Conda: environment management

Detailed methodological references are collected in `docs/methodological_references.md`.
