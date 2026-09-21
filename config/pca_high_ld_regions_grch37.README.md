# PCA high-LD regions (GRCh37)

This file documents the provenance of
`config/pca_high_ld_regions_grch37.bed1`.

## Source

The intervals are the GRCh37/build-37 high-LD-region list distributed in the
public `GrindeLab/PCA` repository with Grinde et al. (2024):

Grinde KE, Browning BL, Reiner AP, Thornton TA, Browning SR.
*Adjusting for principal components can induce collider bias in genome-wide
association studies.* PLoS Genetics. 2024;20(12):e1011242.
doi:10.1371/journal.pgen.1011242.

Upstream source path:
`GrindeLab/PCA/data/highLD/exclude_b37.txt`

The accompanying publication describes the list as regions with high,
long-range, or otherwise unusual LD collected through a literature review.
Build-specific files are distributed for builds 36, 37, and 38.

## Local representation

`pca_high_ld_regions_grch37.bed1` contains the first three columns of the
upstream build-37 file:

`CHROM  START_BP  END_BP`

There is no header because the file is passed directly to PLINK2 with
`--exclude bed1`.  PLINK2's `bed1` modifier denotes 1-based interval input.

The upstream textual source labels are intentionally not used by the command;
the exact interval provenance remains documented here and in
`docs/methodological_references.md`.

## Analytical role

This mask is used only for a predefined PCA sensitivity branch.  The unmasked
MAF+LD-pruned PCA panel is retained in parallel.  No interval is selected from
the observed Cinque Terre PCA results.
