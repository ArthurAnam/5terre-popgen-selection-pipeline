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
doi:10.1093/bioinformatics/btq559. PMID:20926424. PMCID:PMC3025716.
Full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC3025716/

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

### Additional reference for interpretation and limitations

Conomos MP, Reiner AP, Weir BS, Thornton TA. 2016.
**Model-free Estimation of Recent Genetic Relatedness.**
*American Journal of Human Genetics* 98(1):127-148.
doi:10.1016/j.ajhg.2015.11.022.

Role in this project:
- provides an important qualification to a simplistic interpretation of
  KING-Robust values in structured or inbred samples;
- shows that KING-Robust estimates can be negatively biased when at least one
  individual in a pair is inbred, with stronger bias at higher inbreeding;
- shows that ancestry/admixture differences can also bias KING-Robust estimates;
- therefore supports treating negative KING-Robust estimates as properties of
  the estimator/data context rather than as "negative biological relatedness".

### Interpretation notes for manuscript and reviewer responses

The following statements should be kept explicit when KING results are
reported:

- A negative KING-Robust estimate does **not** mean negative biological
  relatedness. The theoretical pedigree kinship coefficient is not interpreted
  as negative; the empirical KING-Robust estimator can nevertheless return
  negative values.
- KING-Robust performs pairwise inference from the two individuals under
  comparison and was designed to reduce confounding by population structure
  relative to estimators that require a single homogeneous set of allele
  frequencies (Manichaikul et al. 2010).
- Negative estimates do not imply "negative biological relatedness". In this
  project they are retained in the numerical outputs and figures because they
  can carry information about pairwise heterogeneity and estimator behaviour.
  This is a project-specific presentation choice, not a universal KING
  requirement: Manichaikul et al. (2010) truncated negative estimates to zero
  in one graphical display, while also showing analytically that unrelated
  individuals from different populations can have a negative KING-Robust
  expectation. Population structure, ancestry differences and inbreeding can
  affect the estimator (Manichaikul et al. 2010; Conomos et al. 2016).
- IBS0 is the proportion of loci at which two individuals carry opposite
  homozygous genotypes and therefore share zero alleles identical-by-state at
  those loci. The kinship-versus-IBS0 plot is a diagnostic recommended by the
  KING documentation for relationship assessment. IBS0 and KING-Robust
  kinship are not statistically independent: opposite-homozygote counts
  contribute directly to the robust kinship estimator, so the descending
  pattern seen in unrelated pairs should not be interpreted as an independent
  biological correlation.
- The conventional KING boundaries used here are >0.354 for duplicate/MZ,
  [0.177,0.354] for first degree, [0.0884,0.177] for second degree, and
  [0.0442,0.0884] for third-degree screening. Third-degree classification is
  treated as exploratory in this project.
- Absence of pairs above 0.0442 should be described as absence of evidence for
  close relatives at the conventional KING third-degree-or-closer threshold,
  not as proof that all individuals are genealogically unrelated at arbitrary
  depth.

### Current Cinque Terre KING result

Using KING 2.3.2 with `--kinship` on the stable pre-HWE 46-sample target
dataset, all 46 choose 2 = 1,035 pairwise comparisons were obtained.

Observed summary:
- maximum KING-Robust kinship estimate: 0.0342;
- minimum estimate: -0.0501;
- no pair reached the exploratory 0.0442 third-degree threshold;
- no duplicate/MZ, first-degree or second-degree candidate pair was detected.

Current QC interpretation:
- no additional sample exclusion is supported by the KING relatedness screen;
- because the retained sample set remains unchanged at 46 individuals, no
  HWE recalculation is required as a consequence of relatedness QC;
- `--related` or `--ibdseg` are not required as routine follow-up in the
  absence of a close or ambiguous candidate pair. They remain available for a
  future IBD-specific scientific question;
- LD pruning is not introduced for KING. This follows the KING documentation,
  which advises retaining good QC-passed SNPs and does not recommend LD pruning
  for relationship inference; the original KING paper likewise reports that
  inference was not impacted by LD structure in the large-sample setting.

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

## PCA / population structure

### Published references reviewed

Patterson N, Price AL, Reich D. 2006.
**Population Structure and Eigenanalysis.**
*PLoS Genetics* 2(12):e190.
doi:10.1371/journal.pgen.0020190. PMCID:PMC1713260.

Role in this project:
- methodological basis for using PCA/eigenanalysis to describe genetic
  population structure;
- supports interpretation of principal components as axes of genome-wide
  genetic variation rather than as predefined discrete population labels;
- notes that the method can be adapted to markers in LD, but does not by itself
  define a universal LD-pruning threshold for every dataset.

Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA, Reich D. 2006.
**Principal components analysis corrects for stratification in genome-wide
association studies.**
*Nature Genetics* 38(8):904-909.
doi:10.1038/ng1847. PMID:16862161.

Role in this project:
- foundational reference for EIGENSTRAT/smartpca-style PCA on genome-wide SNP
  data;
- supports use of principal components to capture ancestry/population
  structure;
- does not justify importing a fixed project-specific MAF or LD threshold
  without checking the current dataset and software documentation.

Sazzini M, Abondio P, Sarno S, et al. 2020.
**Genomic history of the Italian population recapitulates key evolutionary
dynamics of both Continental and Southern Europeans.**
*BMC Biology* 18:51.
doi:10.1186/s12915-020-00778-4. PMID:32438927. PMCID:PMC7243322.

Role in this project:
- especially close methodological precedent because it analyzes high-coverage
  Italian WGS together with external population-reference data;
- used smartpca/EIGENSOFT for population-structure analyses;
- for its PCA-oriented merged dataset, removed one variant from pairs with
  r^2 > 0.2 in windows of 50 SNVs advanced by five SNVs;
- therefore provides a study-specific rationale for adopting
  `--indep-pairwise 50 5 0.2` as the primary PCA LD-pruning strategy here;
- the published PCA dataset was lower density than the original WGS callset, so
  the number and genomic distribution of retained markers must still be audited
  after harmonization in the present project.

### EIGENSOFT / smartpca documentation reviewed

The EIGENSOFT POPGEN documentation records an internal LD-filtering option
(`killr2`) whose default is NO; when enabled, `r2thresh`, `r2genlim` and
`r2physlim` control the filter. The documentation also provides explicit
projection and outlier-handling options.

Project implication:
- LD handling is an explicit PCA-branch decision and is not inherited from KING;
- the primary PCA marker filter is now fixed at MAF >= 0.05 followed by
  PLINK-style LD pruning with 50-SNP windows, 5-SNP steps and r^2 > 0.2,
  following the closely related Sazzini et al. 2020 Italian population study;
- the number and genomic distribution of markers retained after pruning will
  be recorded as an audit rather than assumed to be adequate;
- number of PCs, joint-PCA versus projection strategy, and automatic
  outlier-removal settings remain open until the harmonized Cinque Terre +
  1000 Genomes EUR dataset is inspected.

### 1000 Genomes EUR panel inventory confirmed locally

The local Phase 3 sample metadata file contains 2,504 individuals plus one
header line. Selecting `GROUP == EUR` yields 503 individuals distributed as:

- CEU: 99
- FIN: 99
- GBR: 91
- IBS: 107
- TSI: 107

All 503 EUR sample IDs were confirmed to be present in the chromosome 1 Phase 3
VCF. Before constructing the reproducible EUR subset, the same identity check
will be extended across chromosomes 1-22.

The previously generated local PGEN named
`1KG_EUR.QCcore_mind0.05_geno0.05_alpha0.05.pgen` will not be used as the
primary reference input because it already incorporates historical filtering
choices. The new workflow will instead derive the EUR reference directly from
the original Phase 3 chromosome VCFs and apply the project's documented
harmonization and branch-specific filters explicitly.

## ROH / autozygosity references under review

Howrigan DP, Simonson MA, Keller MC. 2011.
**Detecting autozygosity through runs of homozygosity: a comparison of three
autozygosity detection algorithms.**
*BMC Genomics* 12:460.
doi:10.1186/1471-2164-12-460.

Relevant points for this project:
- PLINK's sliding-window ROH approach was sensitive to marker LD;
- in the study's SNP-array-like data, the authors removed variants with
  MAF < 0.05 before LD pruning and found that LD-pruned data improved
  autozygosity detection;
- these results are important but were developed using SNP-chip-like marker
  densities rather than modern high-density WGS, so they do not yet freeze our
  ROH MAF/pruning settings.

Pemberton TJ, Absher D, Feldman MW, Myers RM, Rosenberg NA, Li JZ. 2012.
**Genomic Patterns of Homozygosity in Worldwide Human Populations.**
*American Journal of Human Genetics* 91(2):275-292.
doi:10.1016/j.ajhg.2012.06.014. PMID:22883143. PMCID:PMC3415543.

Relevant points for this project:
- demonstrates that short, intermediate and long ROH can reflect different
  timescales/processes, from background LD to population history and recent
  parental relatedness;
- used a likelihood-based ROH method on a common intersected SNP panel rather
  than PLINK `--homozyg`, so it informs interpretation but does not directly
  prescribe PLINK-WGS parameters.

Project implication:
- MAF >= 0.05 is now fixed for the PCA branch and already planned for the
  selection/LASSI branch, but it is **not** promoted to a destructive global
  filter on the harmonized master dataset;
- ROH input filtering remains branch-specific and will be finalized after
  reviewing PLINK ROH documentation and WGS-specific marker-density effects.

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
