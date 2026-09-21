# Methodological references and decision log

This file records the literature and software documentation used to justify
analytical choices in the Cinque Terre population-genomics workflow.

The purpose is to keep a direct link between each implemented decision and the
source used to support it. A distinction is made between published literature,
software documentation, and project-specific design decisions.

## KING relatedness

### Published reference

Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM. 2010.
**Robust relationship inference in genome-wide association studies.**
*Bioinformatics* 26(22):2867-2873.

Role in this project:
- primary methodological reference for the KING-Robust kinship estimator;
- supports use of pairwise kinship coefficients for close-relationship inference;
- supports robustness of the estimator to population structure.

### Software documentation

Chen WM. **KING Tutorial: Relationship Inference.** KING documentation,
last updated 28 July 2023.

Role in this project:
- PLINK binary BED/BIM/FAM input;
- use of `--kinship`;
- recommendation not to LD-prune SNPs before KING inference;
- recommendation not to remove otherwise good QC-passed SNPs solely for KING;
- kinship-coefficient boundaries used for duplicate/MZ, first-degree,
  second-degree and exploratory third-degree screening;
- recommendation to inspect kinship versus IBS0;
- interpretation that `--kinship` is strongest for close relationships,
  with third-degree inference treated more cautiously.

## Project-specific KING decisions

The following choices are explicit workflow-design decisions rather than direct
requirements of the KING documentation:

- KING is run on the stable target-QC dataset before final HWE filtering.
  Rationale: relatedness is treated as a sample-level QC/property of the cohort.
  If a close relative were removed, cohort-dependent site statistics, including
  HWE, should be recalculated on the final retained sample set.

- No sample is excluded automatically from KING results. Candidate pairs are
  reviewed first.

- The primary screen is `KING --kinship`. Duplicate/MZ, first-degree and
  second-degree relationships are the main QC targets. The third-degree
  threshold is retained only as an exploratory flag.

- `--related` and `--ibdseg` are not mandatory first-pass analyses. They are
  reserved for follow-up if `--kinship` identifies a close or ambiguous pair
  that would benefit from IBD-segment-based refinement.

## References to add before implementing downstream branches

The following branches must be checked against their original methods papers
and current software documentation before final parameters are frozen:

- PCA / EIGENSOFT: MAF threshold, LD pruning and projection strategy.
- ROH / PLINK: marker-density requirements, LD pruning, MAF handling,
  window/segment parameters, and comparability with 1000 Genomes.
- SHAPEIT phasing: marker filters, genetic map, reference-panel use and
  phasing parameters.
- LASSI: input filtering, MAF handling, phasing requirements, window size and
  shift definition.
- 1000 Genomes reference panel: release, genome build, population definitions
  and harmonization strategy.

For each branch, the final workflow should only promote a working parameter to
a fixed analytical parameter after its supporting documentation has been
reviewed.
