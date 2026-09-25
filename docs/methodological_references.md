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
- the PCA design is now fixed as a joint PCA of the harmonized 46 Cinque Terre
  and 503 1000 Genomes EUR individuals (549 individuals total);
- MAF >= 0.05 is calculated on this joint harmonized dataset before LD pruning;
- LD pruning is then calculated once on the same joint dataset, using 50-SNP
  windows, 5-SNP steps and r^2 > 0.2; there are not separate Cinque Terre and
  1000G pruning lists and no intersection of independently pruned marker sets;
- PLINK2 `--indep-order 1` is specified explicitly so that the pruning order
  follows PLINK 1.x behavior, rather than the newer PLINK2 default order;
- the number and genomic distribution of markers retained after MAF filtering
  and pruning will be recorded as an audit rather than assumed to be adequate;
- number of PCs, explicit long-range-LD-region masking, and automatic
  smartpca outlier-removal settings remain open until this pruned marker set is
  inspected.

#### High/long-range-LD sensitivity design

Price et al. (2008) identified 24 extended-LD regions in European-ancestry
datasets, and the EIGENSOFT documentation recommends excluding the long-range
LD regions from that work before PCA.  Because the Price intervals were
reported as coarse Mb-scale regions in an older genome-build context, this
project does not directly transplant those coordinates onto GRCh37.

Instead, the sensitivity branch uses the explicit GRCh37 catalog distributed
by:

Grinde KE, Browning BL, Reiner AP, Thornton TA, Browning SR. 2024.
**Adjusting for principal components can induce collider bias in genome-wide
association studies.**
*PLoS Genetics* 20(12):e1011242.
doi:10.1371/journal.pgen.1011242. PMID:39680601. PMCID:PMC11684764.

Grinde et al. assembled high, long-range, or otherwise unusual LD regions from
an extensive literature review and distribute coordinate files for genome
builds 36, 37 and 38.  The project's
`config/pca_high_ld_regions_grch37.bed1` is the build-37 interval list from
that resource.

The primary preprocessing already generated without a region mask is retained:
7,548,844 harmonized SNPs -> 4,878,327 SNPs after joint MAF >= 0.05 ->
295,375 SNPs after `--indep-pairwise 50 5 0.2 --indep-order 1`, with all
549 individuals retained.

A second sensitivity branch now starts from the exact same MAF-passing marker
set, excludes the predefined GRCh37 high-LD intervals, and then reruns the
identical LD-pruning parameters.  No region is added or removed in response to
the observed Cinque Terre PCA.  The masked and unmasked PCA results will both
be retained and compared before deciding how the sensitivity analysis should be
presented in the paper; the choice will not be based on which plot appears more
favorable.

The high-LD sensitivity preprocessing completed with 4,878,327 MAF-passing
SNPs. The 18 predefined GRCh37 intervals removed 147,894 SNPs before pruning,
leaving 4,730,433 markers for the masked pruning step. The identical
`--indep-pairwise 50 5 0.2 --indep-order 1` procedure then retained 287,815
masked SNPs, compared with 295,375 SNPs in the unmasked branch. Thus the masked
panel contains 7,560 fewer final PCA markers, while both panels retain the same
549 individuals.

#### smartpca execution settings

Both preprocessed marker panels are carried forward to joint smartpca runs on
the same 549 individuals. Population labels are CT, CEU, FIN, GBR, IBS and TSI.

The smartpca parameterization is explicit:
- `numoutevec: 10`, so the first ten PCs are retained for inspection;
- `numoutlieriter: 0`, disabling automatic iterative sample removal;
- `usenorm: YES`;
- `altnormstyle: NO`, matching EIGENSTRAT-style normalization and the
  EIGENSOFT smartpca example/wrapper behavior;
- `familynames: NO`, preserving the original sample IDs when PACKEDPED/PLINK
  input is read;
- `fastmode: NO`, so the production run uses the exact PCA calculation;
- `numchrom: 22`.

The EIGENSOFT POPGEN documentation states that the default number of output
eigenvectors is ten and that setting `numoutlieriter` to zero disables
outlier removal. Automatic outlier removal is disabled here because the
analysis is explicitly intended to inspect population structure and potential
extreme individuals rather than silently remove them before review.

No `poplistname` is supplied: this is a joint PCA, not a reference-axis
projection. The masked and unmasked smartpca runs differ only in their SNP
panels.

#### Completed joint smartpca run and planned robustness comparison

Both smartpca runs completed with the expected 549 individuals:
CT 46, CEU 99, FIN 99, GBR 91, IBS 107 and TSI 107.

The leading eigenvalues were very similar between preprocessing panels:
- unmasked PC1/PC2/PC3: 4.173998, 2.085995, 1.692304;
- high-LD-masked PC1/PC2/PC3: 4.142286, 2.072798, 1.680901.

The corresponding variance percentages were also similar:
- unmasked: 0.761678%, 0.380656%, 0.308815%;
- high-LD-masked: 0.755892%, 0.378248%, 0.306734%.

These eigenvalue similarities are descriptive only. They do not by themselves
establish that individual coordinates or population placement are unchanged.
The workflow therefore performs a sample-matched quantitative comparison of the
two PCA coordinate sets. For each of the first ten PCs it records the Pearson
correlation, absolute correlation, sign needed for alignment, and aligned RMSE.
It also computes the full 10-by-10 cross-panel correlation matrix to detect
possible component swaps or rotations, population-centroid shifts for PC1-PC3,
and CT-to-reference centroid distances. PC signs are explicitly aligned before
coordinate differences are interpreted because eigenvector sign is arbitrary.

The masked-versus-unmasked interpretation will be based on these quantitative
comparisons and the corresponding plots, rather than on visual preference for
one PCA figure.

The coordinate comparison has now completed. For PC1-PC3, absolute same-PC
correlations were 0.999725, 0.998224 and 0.997074, respectively, with mean
absolute correlation 0.998341. The first three components therefore show very
high coordinate concordance between the masked and unmasked marker panels.
The maximum population-centroid shift across PC1-PC3 was 0.005218.

CT-to-reference centroid distances across PC1-PC3 were also nearly unchanged:
CEU 0.157634 -> 0.157442; FIN 0.175689 -> 0.175663; GBR 0.161689 -> 0.162084;
IBS 0.160217 -> 0.160247; TSI 0.160447 -> 0.160297. These differences are
small relative to the corresponding centroid distances.

The higher PCs are less stable one-to-one. PC4-PC6 remain strongly correlated
with their same-index counterparts, but PCs 7-10 show lower same-index
correlations and cross-match to neighboring PCs. This pattern is consistent
with rotation/reordering among later components whose eigenvalues are close;
it is not interpreted as a contradiction of the PC1-PC3 robustness result.

For visual review, the comparison workflow also generates side-by-side population plots
for PC1-PC2 and PC2-PC3. The masked coordinates are sign-aligned to the
unmasked solution and both panels use identical axis limits, so apparent
differences are not introduced by arbitrary eigenvector sign changes or
automatic rescaling.

To limit redundant disk use, PCA figures are stored in a single raster format (PNG at 300 dpi) rather than duplicated as both PNG and PDF. Vector output can be regenerated later at manuscript-submission stage without rerunning smartpca.

Paper presentation remains pending visual review of the PC1-PC2 and PC2-PC3
plots. Because the primary structure is quantitatively robust to the mask, the
choice of which panel is shown as the main figure can be based on a priori
methodological framing rather than on which version gives a preferred result.


### Exploratory population grouping diagnostic

A lightweight exploratory branch tests whether the current population labels
(CT, CEU, FIN, GBR, IBS, TSI) show clear genome-wide similarity patterns.
It is an explicit target and does not alter any population label elsewhere in
the workflow.

For speed and consistency, the diagnostic reuses the existing joint
MAF >= 0.05, predefined high-LD-masked and LD-pruned PCA marker panel
(287,815 SNPs; 549 individuals). Pairwise differentiation is estimated with
PLINK2 using the Hudson FST estimator. PLINK2 documents Hudson as the default
FST method and cites Bhatia et al. 2013 for this estimator.

The pairwise FST values are reported directly and also used as distances for an
average-linkage hierarchical dendrogram. Negative finite FST estimates, if any,
are retained in the raw outputs but truncated to zero only when constructing
the clustering distance, because a negative sampling estimate is not a
biological negative distance.

This branch is deliberately descriptive: it ranks the closest population pairs
and visualizes their hierarchical similarity, but it does not define a
threshold for declaring populations identical and does not automatically merge
groups. If a proposed pooling becomes important to downstream inference, a
separate uncertainty/stability analysis should be added before formal use.

Observed exploratory FST results
--------------------------------

The completed quick diagnostic produced the following pairwise Hudson FST
ranking (lowest to highest): CEU-GBR 0.00031653; IBS-TSI 0.0015252;
CEU-IBS 0.00230597; GBR-IBS 0.0025142; CEU-TSI 0.00346362;
GBR-TSI 0.00385377; CT-TSI 0.0053019; CT-IBS 0.00603848;
CEU-FIN 0.00618549; FIN-GBR 0.00659861; CEU-CT 0.0077535;
CT-GBR 0.00813886; FIN-IBS 0.010168; FIN-TSI 0.0116828; and
CT-FIN 0.0158024.

Average-linkage clustering first joined CEU with GBR (height 0.00031653),
then IBS with TSI (0.0015252), then those two reference clusters
(0.00303439). CT joined that four-population cluster at 0.006808185, and FIN
joined last at 0.01008746. This is retained as a descriptive exploratory
result only; it is not used to merge labels in downstream analyses.

Primary estimator reference:
Bhatia G, Patterson N, Sankararaman S, Price AL. 2013. Estimating and
interpreting FST: The impact of rare variants. Genome Research 23:1514-1521.
doi:10.1101/gr.154831.113.

### 1000 Genomes EUR panel inventory confirmed locally

The local Phase 3 sample metadata file contains 2,504 individuals plus one
header line. Selecting `GROUP == EUR` yields 503 individuals distributed as:

- CEU: 99
- FIN: 99
- GBR: 91
- IBS: 107
- TSI: 107

All 503 EUR sample IDs were confirmed to be present in every autosomal Phase 3
VCF from chromosomes 1-22.

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


### 1000 Genomes EUR downstream QC design

The EUR reference branch is derived from the original 1000 Genomes Phase 3
chromosome VCFs rather than from previously filtered local derivative files.

Confirmed EUR sample inventory:
- CEU: 99
- FIN: 99
- GBR: 91
- IBS: 107
- TSI: 107
- total EUR samples: 503

All 503 EUR sample IDs were verified to be present in the local Phase 3
chromosome VCF collection. All autosomal VCFs are indexed.

Downstream QC principles for the EUR subset:
- retain autosomes only;
- retain biallelic SNPs;
- recalculate cohort-dependent AC, AN and missingness after restricting to the
  503 EUR individuals;
- remove variants that become monomorphic in the EUR subset;
- remove sites with missingness > 5%;
- evaluate individual missingness > 5% and repeat cohort-dependent site QC only
  if an individual is removed;
- do not use pooled-EUR HWE as a filtering criterion;
- HWE may be inspected descriptively within CEU, FIN, GBR, IBS and TSI
  separately, without removing variants on that basis;
- strand-ambiguous SNPs are handled during harmonization with the Cinque Terre
  dataset.

Genome-wide EUR-subset QC audit across chromosomes 1-22:
- biallelic SNPs inspected: 77,818,346
- variants monomorphic after restriction to EUR: 55,685,456
- sites with any missing genotype: 0
- sites with missingness > 5%: 0
- polymorphic biallelic SNPs before strand-ambiguity removal: 22,132,890
- palindromic A/T or C/G SNPs removed: 3,381,797
- polymorphic non-palindromic SNPs retained for harmonization: 18,751,093

Because no missing genotype was observed at any inspected biallelic SNP, the
site-missingness filter removes zero variants and individual missingness is
necessarily zero in this EUR subset. The main reductions at this stage are
therefore caused by variants becoming monomorphic after restriction from the
global Phase 3 cohort to the 503 EUR individuals and by the planned removal of
strand-ambiguous palindromic SNPs.


### ROH marker-density audit before parameter freezing

The ROH branch now fixes only the MAF threshold (MAF >= 0.05) and deliberately
does not yet freeze LD pruning, minimum SNP count, minimum physical length,
maximum gap, or density settings. The first reproducible ROH step is a marker
density audit on the exact harmonized 549-sample CT + 1000G EUR dataset.

For this initial audit, MAF >= 0.05 is calculated jointly across the 549
harmonized samples, matching the already observed joint common-variant count of
4,878,327. No LD pruning is applied at this stage. The audit reports, by
chromosome, the number of retained markers, physical span, marker density,
inter-marker gap distribution, and the fractions of consecutive-marker gaps
exceeding 50, 100, 250, 500 and 1000 kb.

This is intentionally a preprocessing diagnostic rather than a ROH call.
PLINK 1.9 documents that ROH results depend on minimum SNP count, minimum
physical span, mean SNP density and maximum internal gap, while Howrigan et al.
(2011) showed that LD treatment and marker density can materially affect PLINK
ROH detection. Therefore these parameters will be chosen only after the
observed density of the present WGS-derived common-marker panel is known.

The joint MAF scope used for this density audit is provisional for the final
cross-population ROH comparison. Before the production ROH call, the project
will explicitly assess whether the common comparison panel should require the
MAF threshold jointly or within each population, so population-specific allele
frequency differences are not silently converted into a marker-ascertainment
difference.

#### Completed ROH marker-density audit

The audit completed on the expected 549-sample harmonized dataset and retained
exactly 4,878,327 variants after joint MAF >= 0.05. Across the terminal-marker
spans of the 22 autosomes, the panel covers 2,793,824,144 bp and contains
1,746.11 markers/Mb, corresponding to a genome-wide mean inter-marker gap of
572.70 bp.

The panel is therefore far denser than a conventional SNP-array panel. The
chromosome-specific median gap ranges from 247 to 318 bp and the chromosome-
specific 99th percentile is approximately 3.1-4.0 kb. Large gaps are rare:
355 of 4,878,305 consecutive-marker gaps exceed 50 kb, 201 exceed 100 kb,
61 exceed 250 kb, 35 exceed 500 kb, and 21 exceed 1 Mb. The largest observed
gap is 21,133,173 bp.

Interpretation:
- minimum SNP-count thresholds cannot be interpreted independently of physical
  length on this WGS-derived panel, because even a short physical interval can
  contain many common markers;
- the final ROH definition should therefore include an explicit physical-length
  threshold and should not rely on an SNP-count threshold alone;
- the handful of very large gaps are localized coverage/assembly exceptions and
  should be prevented from bridging otherwise separate ROH by an explicit
  maximum-gap criterion;
- LD pruning remains the main unresolved preprocessing choice. Howrigan et al.
  supports light-to-moderate pruning in SNP-array-like data, while later work
  shows that pruning effects can be population- and dataset-dependent. The next
  workflow step should therefore compare an unpruned common-marker call with a
  light/moderate LD-pruned sensitivity analysis before freezing the production
  ROH settings.

### ROH decision framework after literature review

The ROH literature review is recorded here rather than only in generated
README files. Generated README files explain individual result directories;
`config/config.yaml` records executable parameter choices; this document
records the biological and methodological rationale behind those choices; and
`docs/analysis_records/` stores compact observed numerical results that would
otherwise live only under ignored `results/` paths.

The workflow deliberately uses analysis-specific LD treatment. The PCA branch
uses strong pairwise-r2 pruning because PCA is intended to summarize broad
population structure without allowing local blocks of highly correlated
markers to dominate eigenvectors. That PCA marker panel is not automatically
appropriate for ROH detection, where extended homozygosity and local LD are
part of the biological/statistical problem being measured.

For ROH, two literature-supported approaches are therefore kept distinct:

1. Primary candidate: the common harmonized MAF >= 0.05 marker panel without
   LD pruning, with physical ROH length intended to be the principal protection
   against short background-LD tracts. This follows the population-history /
   FROH tradition in which long ROH are analyzed as genomic segments rather
   than forcing approximate marker independence.
2. Methodological sensitivity: PLINK 1.9 VIF pruning using
   `--indep 50 5 10`. This reproduces the light-pruning procedure evaluated
   by Howrigan et al. (2011) rather than replacing it with
   `--indep-pairwise`. The VIF threshold of 10 corresponds to multiple
   R-squared = 0.90 in the single-predictor special case, but the algorithms
   are not equivalent.

The existing PCA pruning (`--indep-pairwise 50 5 0.2`) is intentionally not
reused for ROH. Reusing it merely for pipeline uniformity would impose a much
stronger, differently defined LD filter than the ROH-specific methodological
precedent and could remove informative homozygous sequence.

A minimum physical ROH length of 1.5 Mb is retained as a strong candidate, not
yet a frozen production parameter. McQuillan et al. (2008) showed that ROH
>=1.5 Mb distinguish European populations with different isolation histories.
The northeastern-Italian isolate study also defined gROH using a 1.5 Mb
minimum, and later human ROH work commonly treats this scale as a useful
boundary between shorter background-homozygosity tracts and longer segments
more informative about autozygosity/population history. Because the present
unpruned panel contains ~1,746 common markers/Mb, physical length is more
informative than importing a historical SNP-count threshold from a much
sparser array.

The light-VIF density audit is now complete. PLINK 1.9
`--indep 50 5 10` retained 814,440 of the 4,878,327 common MAF>=0.05
markers. Across the autosomal terminal-marker spans this corresponds to
291.52 markers/Mb and a mean inter-marker gap of 3.43 kb. Only 207 consecutive
marker gaps exceed 100 kb, 35 exceed 500 kb, and 21 exceed 1 Mb.

Thus, even after light VIF pruning, the panel remains dense: a 1.5 Mb interval
contains roughly 437 retained markers on average. A historical 50-65 SNP
minimum is therefore a safety floor in this dataset, not the biological
definition of a long ROH. Physical length remains explicit.

The next step is a prespecified two-framework sensitivity analysis. The
population-history candidate uses the unpruned common panel and calls ROH
directly at >=1.5 Mb with PLINK settings 50 SNP/window, minimum 50 SNP per called ROH, density
50 kb/SNP, gap 1000 kb, 5 missing/window, 1 heterozygote/window, and window
threshold 0.05. The light-VIF sensitivity uses `--indep 50 5 10`, then the
Howrigan/UK-Biobank-style 65 SNP/window, minimum 65 SNP per called ROH, density 200 kb/SNP,
gap 500 kb, 3 missing/window, 0 heterozygotes/window, threshold 0.05, and a
minimal 10-kb calling floor. For comparability, only segments >=1.5 Mb from
that sensitivity call are retained in downstream summaries.

This design intentionally does not compute FROH yet. First, N_ROH, summed ROH
length, mean/median ROH length and maximum ROH length are compared at the same
>=1.5 Mb biological scale. The FROH denominator will be frozen and documented
separately after the production ROH definition is selected.

#### Completed long-ROH call sensitivity

The two prespecified frameworks completed successfully on all 549 individuals.
The unpruned population-history candidate produced a mean summed ROH burden of
37.55 Mb in CT, compared with 20.60 Mb in FIN, 14.42 Mb in IBS, 12.91 Mb in
GBR, 9.43 Mb in CEU and 8.89 Mb in TSI. Under the light-VIF/Howrigan-derived
sensitivity, the corresponding means were 28.19 Mb in CT, 7.05 Mb in FIN,
4.79 Mb in IBS, 4.00 Mb in GBR, 2.44 Mb in TSI and 2.37 Mb in CEU.

Thus CT has the highest mean burden of ROH >=1.5 Mb under both definitions.
Absolute burden decreases under light VIF pruning in every population, but the
reduction is much smaller for CT: the light-VIF framework retains ~75.1% of
the CT mean total burden versus ~25.1-34.2% in the reference populations.
This is treated as robustness evidence, not as a standalone demographic proof,
because the two frameworks differ jointly in LD treatment and several linked
PLINK parameters.

At the individual level within CT, agreement is especially strong. Total ROH
burden has Pearson r=0.993 and Spearman rho=0.977 between frameworks; N_ROH
has Pearson r=0.909 and Spearman rho=0.911. Hence individuals with high versus
low long-ROH burden remain ordered very similarly even though absolute segment
burden changes.

The next sensitivity deliberately returns to the unpruned primary candidate
and holds marker panel, MAF, minimum physical length, SNP count, density,
missingness and window threshold fixed. Only maximum internal gap
(100/500/1000 kb) and tolerated heterozygotes per window (0/1) are varied in a
3x2 design. This isolates the two remaining technical choices before the
production ROH definition and FROH denominator are frozen.

#### Completed gap x heterozygote sensitivity

The 3x2 sensitivity completed successfully. With one heterozygote allowed per
50-SNP window, changing the maximum gap from 1000 to 500 kb had only a small
effect. In CT, mean total ROH burden changed from 37.55 to 36.58 Mb (97.4%
retained), with Pearson r=0.9995 and Spearman rho=0.9975 at the individual
level. Maximum ROH length in CT was perfectly correlated between these two
gap settings. Gap=500 kb is therefore selected as the primary candidate: it
is more conservative against bridging across long marker-free intervals than
1 Mb while preserving essentially the same long-ROH signal.

In contrast, setting the tolerated heterozygote count to zero had a large
effect. At gap=500 kb, mean total burden with zero versus one heterozygote per
window retained ~60.4% in CT but only ~20.5% in CEU, 29.6% in FIN, 26.9% in
GBR, 28.9% in IBS and 23.5% in TSI. This population/dataset asymmetry makes a
zero-heterozygote production rule undesirable for the present mixed-source
CT+1000G comparison. Human WGS work has emphasized that isolated heterozygous
calls can split true ROH because of sequencing/genotyping error, and PLINK's
heterozygote tolerance exists to mitigate this problem. Allowing one
heterozygote per window is therefore selected as the primary candidate. This
choice is also conservative with respect to a CT-vs-reference contrast because
it increases the retained ROH burden proportionally more in the reference
populations than in CT.

The resulting technical production candidate is: unpruned common panel,
ROH >=1.5 Mb, 50 SNP/window, minimum 50 SNP per called ROH, density <=50 kb/SNP, gap <=500 kb,
5 missing genotypes/window, 1 heterozygote/window and window threshold 0.05.

One marker-ascertainment question remains before this definition is frozen:
the current common panel uses MAF>=0.05 calculated jointly across all 549
samples. A final audit compares this with the intersection of SNPs having
MAF>=0.05 within each of CT, CEU, FIN, GBR, IBS and TSI.

#### Completed MAF-scope density audit

The strict all-population intersection retains 4,142,839 of 4,878,327 joint
MAF>=0.05 variants (84.92%). Considered separately, each population retains
94.10-96.73% of the joint panel: CT 94.10%, FIN 94.17%, GBR 95.22%, TSI
96.29%, IBS 96.53% and CEU 96.73%. The larger 15.08% loss in the intersection
therefore reflects the accumulation of population-specific frequency
differences rather than one population losing a very large fraction alone.

The strict intersection remains extremely dense: 1,482.88 markers/Mb with a
mean inter-marker gap of 674 bp. Large-gap structure is nearly unchanged
relative to the joint panel: 208 gaps exceed 100 kb, 37 exceed 500 kb and 22
exceed 1 Mb, versus 201, 35 and 21 respectively in the joint-MAF panel.

The direct MAF-scope call sensitivity is now complete. With otherwise
identical ROH parameters, the strict all-population MAF intersection produced
slightly higher absolute long-ROH burden than the joint-MAF panel (CT mean
37.65 versus 36.58 Mb; FIN 20.81 versus 19.41 Mb; IBS 14.85 versus 13.63 Mb;
GBR 13.44 versus 11.94 Mb; CEU 9.65 versus 8.62 Mb; TSI 9.16 versus 8.14 Mb).

Within CT, total long-ROH burden has Pearson r=0.9990 and Spearman rho=0.9925
between marker ascertainment schemes; N_ROH has Pearson r=0.9903 and Spearman
rho=0.9810, and maximum ROH length has Pearson r=0.9909 and Spearman rho=0.9965.
Joint MAF>=0.05 is therefore frozen as the primary marker panel because it is
less restrictive and does not require every differentiated SNP to be common in
every reference population. The all-population MAF>=0.05 intersection remains
a prespecified sensitivity analysis.

#### Final production ROH definition and FROH

The final production ROH definition was frozen on 25 September 2026 after
dataset-specific sensitivity analyses of missing-genotype tolerance,
heterozygote tolerance, maximum inter-marker gap, MAF ascertainment and
LD pruning.

Earlier ROH candidate settings documented above are retained as methodological
provenance but are superseded by this final calibration.

The production marker panel contains 549 individuals (46 Cinque Terre and
503 1000 Genomes Phase 3 EUR individuals) and 4,878,327 autosomal SNPs after
joint MAF >=0.05 filtering. The same marker panel is used for all six
populations. No LD pruning is applied for the primary ROH call.

The final PLINK 1.9 `--homozyg` settings are:

- `--homozyg-window-snp 50`
- `--homozyg-snp 50`
- `--homozyg-kb 1500`
- `--homozyg-density 50`
- `--homozyg-gap 250`
- `--homozyg-window-missing 1`
- `--homozyg-window-het 2`
- `--homozyg-window-threshold 0.05`

No global `--homozyg-het` limit is imposed.

The distinction between the two SNP-count parameters is explicit.
`--homozyg-window-snp 50` specifies the local sliding-window size, whereas
`--homozyg-snp 50` is a minimum SNP count for the final called segment.
The latter is non-binding in the final dense WGS-derived panel: under the
final gap=250 kb configuration the smallest observed ROH contains well above
50 SNPs. A Lencz/Purfield-style minimum-marker calculation gave a diagnostic
value of approximately 64 SNPs, but this value was not imposed as an
additional threshold because the calculation assumes properties that do not
fully represent the LD structure of the present WGS panel. The empirical final
callset demonstrates that the 50-SNP technical floor does not determine the
called segments.

A minimum physical length of 1.5 Mb remains the primary biological ROH
definition. ROH >=5 Mb are additionally summarized as a secondary
length-specific endpoint because the sensitivity analyses showed that the
long-ROH burden is especially stable and is informative about the long-run
component of autozygosity.

##### Missing-genotype tolerance

The production setting permits at most one missing genotype per 50-SNP
sliding window.

A direct sensitivity comparison of 0, 1 and 5 missing calls per window showed
that the effect was confined primarily to the Cinque Terre samples because the
retained 1000 Genomes EUR panel has no missing genotypes at the harmonized
sites, whereas CT retains a very small amount of missingness.

Allowing zero missing calls fragmented some CT ROH. Allowing one versus five
missing calls produced only a small change in CT total burden. The final
choice of one missing call per window therefore provides limited tolerance for
the observed CT missingness without retaining the substantially more permissive
PLINK default of five.

##### Heterozygote tolerance

The production setting permits at most two heterozygous calls per 50-SNP
sliding window.

Sensitivity analyses from zero through four tolerated heterozygotes showed a
large increase from zero to one and from one to two, followed by progressively
smaller gains in long-ROH burden. Detailed inspection demonstrated that
`window-het=1` can fragment extremely long, overwhelmingly homozygous regions.

For example, in NA20585 on chromosome 8, the setting with one tolerated
heterozygote split a region into several adjacent ROH, whereas
`window-het=2` reconstructed a single approximately 37.7-Mb tract containing
61,804 SNPs and only 41 heterozygotes overall
(PHET approximately 0.00066). Similar behavior was observed in other audited
long segments.

In contrast, the higher-PHET segments introduced under `window-het=2` in CT
were concentrated among short ROH close to the 1.5-Mb lower boundary.
Among CT segments with PHET >=0.02 under this setting, none reached 5 Mb.
The >=5-Mb burden therefore showed substantially greater stability than the
short-ROH component.

A global `--homozyg-het` limit is intentionally not used. A fixed absolute
number of heterozygotes across the complete ROH would become increasingly
stringent with segment length and could reject very long segments that remain
overwhelmingly homozygous. Local heterozygote control is instead imposed by
the sliding-window criterion.

##### Maximum mean SNP density

The production setting retains `--homozyg-density 50`, corresponding to a
maximum mean spacing of 50 kb per SNP within a called ROH.

This threshold is empirically non-binding in the final WGS-derived callset.
Among the 3,933 production ROH, the median observed DENSITY value is
0.658 kb/SNP, the 95th percentile is 1.4774 kb/SNP, the 99th percentile is
5.567 kb/SNP and the maximum is 11.383 kb/SNP.

A direct sensitivity analysis reran the otherwise identical final production
call with `--homozyg-density 100` and `--homozyg-density 1000`. Both runs
returned exactly 3,933 ROH, and the complete segment identities
(FID, IID, chromosome, start, end and SNP count) were identical to the
density=50 production call.

Therefore `--homozyg-density 50` is retained as a technical guardrail rather
than a parameter determining segment inclusion in the present dataset. The
observed marker density of all final ROH lies comfortably inside this limit.

##### Maximum inter-marker gap

The final maximum gap is 250 kb.

A fresh sensitivity analysis compared 100, 250 and 500 kb while holding the
final missingness and heterozygote settings fixed. In CT, mean total ROH burden
was approximately 39.39, 40.37 and 41.42 Mb for gap thresholds of 100, 250 and
500 kb, respectively. Mean burden from ROH >=5 Mb was approximately 18.23,
18.53 and 18.66 Mb.

Thus 250 kb retains almost all of the long-ROH burden recovered at 500 kb,
while avoiding frequent bridging across 250-500-kb marker-sparse intervals.
The 500-kb setting produced hundreds of segments containing an internal gap
above 250 kb, including recurrent calls in known marker-sparse/pericentromeric
regions. Conversely, 100 kb was sufficiently restrictive to fragment some
otherwise highly homozygous long tracts.

The 250-kb value is therefore an empirical compromise between fragmentation
and bridging rather than an imported generic default.

##### MAF-scope sensitivity

Joint MAF >=0.05 is used for the primary cross-population marker panel.

The final-parameter sensitivity compared three marker definitions:

1. joint MAF >=0.05 across all 549 individuals;
2. the intersection of SNPs having MAF >=0.05 separately in CT, CEU, FIN,
   GBR, IBS and TSI;
3. the harmonized panel without a MAF filter.

The primary joint-MAF panel contains 4,878,327 SNPs. The strict six-population
intersection contains 4,142,839 SNPs, while the no-MAF harmonized panel
contains 7,548,844 SNPs.

Within CT, total ROH burden was highly concordant with the joint-MAF result:
Pearson r=0.9986 for the strict all-population intersection and r=0.9982 for
the no-MAF panel. For burden contributed by ROH >=5 Mb, the corresponding
Pearson correlations were 0.9980 and 0.9943.

Mean CT burden from ROH >=5 Mb was 18.53 Mb with joint MAF >=0.05,
18.83 Mb with the strict all-population intersection and 17.54 Mb without a
MAF filter.

These results show that the main long-ROH pattern is not created by the chosen
MAF ascertainment. Joint MAF >=0.05 is retained as primary because it provides
one common marker panel for all individuals without requiring every
differentiated locus to be common independently in all six populations.

##### LD-pruning sensitivity

No LD pruning is used in the primary ROH analysis.

A final sensitivity applied the VIF-based procedures considered by
Howrigan et al. (2011): light pruning with `--indep 50 5 10` and moderate
pruning with `--indep 50 5 2`.

The unpruned joint-MAF panel contains 4,878,327 SNPs. Light VIF pruning retains
814,440 SNPs and moderate VIF pruning retains 426,484 SNPs.

With a fixed 50-SNP PLINK sliding window, this thinning changes the physical
scale of the scan substantially. The median physical span of a consecutive
50-SNP window is approximately:

- 23.2 kb in the unpruned panel;
- 144.3 kb after light VIF pruning;
- 282.7 kb after moderate VIF pruning.

The corresponding 95th percentiles are approximately 54.3, 321.1 and
560.6 kb.

Thus LD pruning in this WGS-density dataset is not merely removal of correlated
markers; with fixed marker-count windows it also changes the physical scale of
the ROH algorithm.

Nevertheless, individual CT burden remains highly correlated with the unpruned
analysis. For total ROH burden, Pearson r is 0.9978 under light pruning and
0.9958 under moderate pruning. For ROH >=5 Mb, Pearson r is 0.9878 and 0.9872,
respectively.

Absolute long-ROH burden increases rather than decreases monotonically with
pruning: mean CT burden from ROH >=5 Mb is 18.53 Mb unpruned, 20.92 Mb after
light pruning and 22.57 Mb after moderate pruning. This demonstrates that
marker thinning alters segment definition and is not intrinsically a more
conservative ROH procedure in this dataset.

The unpruned WGS-derived panel is therefore retained for production, while the
VIF-pruned calls serve as methodological sensitivity analyses.

##### FROH definition

Primary FROH is calculated as the summed physical length of all autosomal
ROH >=1.5 Mb divided by 2.77e9 bp.

The same denominator is used for every individual and population. ROH lengths
are retained continuously, and an additional FROH based only on ROH >=5 Mb is
reported as a secondary length-specific endpoint.

#### Final production ROH and FROH results

The final production workflow completed successfully on all 549 individuals
and produced 3,933 ROH >=1.5 Mb.

For CT (n=46):

- mean N_ROH >=1.5 Mb = 10.87;
- median N_ROH >=1.5 Mb = 9;
- mean total ROH burden >=1.5 Mb = 40.37 Mb;
- median total burden >=1.5 Mb = 25.74 Mb;
- mean FROH >=1.5 Mb = 0.01457;
- median FROH >=1.5 Mb = 0.00929;
- mean number of ROH >=5 Mb = 2.00;
- median number of ROH >=5 Mb = 1;
- mean burden from ROH >=5 Mb = 18.53 Mb;
- median burden from ROH >=5 Mb = 6.30 Mb;
- mean FROH >=5 Mb = 0.00669;
- median FROH >=5 Mb = 0.00228.

Mean total ROH burden >=1.5 Mb in the reference populations is 13.05 Mb in
CEU, 24.51 Mb in FIN, 16.92 Mb in GBR, 18.13 Mb in IBS and 12.13 Mb in TSI.

The distribution analysis uses two-sided Mann-Whitney tests for each
prespecified CT-versus-reference comparison, Holm correction across the five
comparisons separately within each endpoint, and Cliff's delta as an effect
size.

For primary FROH >=1.5 Mb, CT differs from CEU
(Holm-adjusted p=2.88e-8; Cliff's delta=0.598), GBR
(p=2.90e-5; delta=0.461), IBS (p=2.90e-5; delta=0.452) and TSI
(p=6.06e-10; delta=0.657). CT versus FIN is not significant after Holm
correction (p=0.215; delta=0.128).

For FROH contributed specifically by ROH >=5 Mb, CT differs from all five
reference populations, including FIN. The Holm-adjusted p-values are
2.60e-14 for CEU, 2.27e-5 for FIN, 1.34e-8 for GBR, 2.11e-6 for IBS and
9.40e-9 for TSI.

The distinction between the two endpoints is biologically relevant to the
interpretation: the total >=1.5-Mb FROH distribution does not show a clear
uniform CT-versus-FIN shift, whereas the very-long-ROH component is elevated
in CT relative to FIN as well as to the other EUR reference populations.

These results support interpretation of the CT signal as having a particularly
strong long-ROH component. ROH length is an indirect marker of demographic
history, however, and the analysis is not used to infer a precise number of
generations or to attribute the pattern to a single historical event.


### Howrigan et al. 2011: specific ROH parameter implications

Howrigan et al. (2011) is a key parameter-tuning reference for PLINK ROH
calling. Important details for this project:

- before ROH calling, they removed SNPs with MAF < 0.05;
- they used PLINK `--indep` (VIF-based pruning), not pairwise-r2 pruning;
- their "light" LD pruning used a 50-SNP window with VIF > 10
  (approximately r2 > 0.9 for a simple pair);
- their "moderate" LD pruning used a 50-SNP window with VIF > 2
  (approximately r2 > 0.5 for a simple pair);
- stronger pruning (e.g. VIF ~1.33, approximately r2 > 0.25) performed worse
  than light-to-moderate pruning in their simulations;
- the best-performing analyses did not improve when one heterozygous genotype
  was allowed inside called ROH; the authors therefore recommended allowing
  zero heterozygotes in called ROH;
- the optimal minimum SNP count depended on the target autozygosity timescale:
  with moderate pruning, ~45-50 SNPs performed best for autozygosity within
  ~20 generations, while ~35 SNPs performed best for ~50 generations;
- with light pruning, ~65 SNPs was preferred for the ~50-generation scenario;
- the corresponding best-performing minimum physical span was roughly
  750 kb in the examples highlighted by the authors.

Caution:
- the study simulated sequence data but then subsampled common variants to
  mimic dense SNP-array data, so these settings are strong guidance rather
  than an automatic final prescription for present-day high-density WGS;
- because the paper used `--indep`, translating its pruning settings into
  `--indep-pairwise` would not be a methodologically exact reproduction.
  If this paper is used as the main ROH-pruning precedent, a VIF-based
  implementation should be considered directly.


### 1000 Genomes EUR harmonization audit

A chromosome-by-chromosome exact-key audit was performed using
`CHR:POS:REF:ALT` after EUR subsetting, restriction to polymorphic biallelic
SNPs, the missingness checks defined above, and removal of palindromic A/T and
C/G sites.

Genome-wide counts:
- Cinque Terre QC SNPs: 9,000,246
- 1000G EUR non-palindromic QC records: 18,751,093
- 1000G EUR unique `CHR:POS:REF:ALT` keys: 18,751,091
- sites shared by chromosome and position: 7,551,493
- exact `CHR:POS:REF:ALT` matches: 7,548,844
- same-position allele mismatches: 2,649

The two duplicated 1000G EUR exact keys were:
- `12:8400000:T:G`
- `17:1144632:C:T`

For each duplicated key, the two 1000G records had identical EUR genotype
vectors (verified by matching SHA256 hashes). Neither locus was present in the
Cinque Terre QC VCF, so these duplicate records do not affect the exact
harmonized intersection.

Harmonization policy:
retain only exact `CHR:POS:REF:ALT` matches and exclude same-position
allele-discordant sites; no strand flipping or allele-rescue procedure is used.


#### Reproducible implementation note

The production harmonization branch derives the 503-person EUR sample list
directly from the local `1000GP_Phase3.sample` metadata file rather than
depending on the temporary manual list used during the audit. Each chromosome
is processed independently to limit disk use. The reference chromosome is
subset to EUR samples, restricted to polymorphic biallelic SNPs, checked for
the verified zero-missingness condition after site QC, and stripped of
palindromic A/T and C/G sites. Exact matches with the final Cinque Terre QC VCF
are then identified with allele-exact `bcftools isec` logic and merged across
samples. Same-position allele-discordant records are not rescued or flipped.
Per-chromosome BCFs and the concatenated BCF are temporary workflow
intermediates; the persistent harmonized master dataset is stored in PLINK2
PGEN/PVAR/PSAM format for reuse by downstream branches.


#### Reproducible run confirmation

The production harmonization workflow reproduced the manual genome-wide audit
exactly and completed with 46 Cinque Terre samples plus 503 1000G EUR samples
(549 total). The persistent harmonized PGEN contains 7,548,844 variants.

The workflow reports 18,751,093 EUR QC records before intersection. This is a
record count, whereas the manual exact-key audit contained 18,751,091 unique
`CHR:POS:REF:ALT` keys because the source callset contains the two duplicated
keys documented above. Neither duplicate locus occurs in the Cinque Terre QC
VCF, so the record-versus-unique-key distinction does not change the final
7,548,844-site harmonized intersection.


### Selection branch: LD-decay calibration for LASSI

The selection scan is now entering a separate branch from the completed
population-structure analyses. The primary LASSI target is the final QCed
Cinque Terre cohort (46 individuals); reference populations are not pooled into
the CT LD calculation.

The original LASSI empirical human protocol based its SNP-delimited window on
the physical interval over which pairwise LD, measured as r^2, decayed below
one third of its value for SNP pairs separated by 1 kb. Harris and DeGiorgio
used 117-SNP windows with a 12-SNP step in their human application; those
numbers are not copied directly because CT marker density and LD structure are
dataset-specific.

For CT, physical LD decay is therefore estimated first. The primary audit uses
MAF >=0.05, matching the working selection-preprocessing threshold already
recorded in the project, while MAF >=0.01 is run as a prespecified sensitivity.
The goal is not to use MAF filtering as a biological definition of selection;
it is to establish whether the physical LD-decay estimate is robust to the
common-variant threshold in this 46-sample cohort.

An exhaustive all-pairs calculation is not appropriate for the WGS-density
panel. Instead, 100 anchor SNPs per autosome are selected deterministically at
approximately even physical positions after excluding the first and last
500 kb of each chromosome's observed marker span. Each anchor is compared
with every eligible SNP within 500 kb using PLINK 1.9 unphased hard-call r^2.
This produces thousands of observations per 1-kb distance bin while keeping
the intermediate pair table tractable. Anchor selection is deterministic, so
the estimate is exactly reproducible.

The baseline LD value is the mean r^2 of pairs separated by 0.5-1.5 kb.
The direct LASSI-style decay estimate is the first subsequent 1-kb bin whose
mean r^2 is below one third of that baseline. A requirement of five
consecutive bins below the threshold is also reported as a stability
diagnostic, but is not substituted silently for the literature-based first
crossing.

This stage deliberately stops at a physical distance. The final LASSI
window size in SNPs will be derived only after the MAF-sensitivity result is
reviewed and the marker panel used for the scan is frozen. The step size will
then be approximately 10% of the final SNP-window size, matching the design of
the original LASSI human scan.

Reference correction: the primary LASSI method is Harris AM & DeGiorgio M
(2020), "A likelihood approach for uncovering selective sweep signatures from
haplotype data", Molecular Biology and Evolution 37:3023-3046,
doi:10.1093/molbev/msaa115. The previously listed Genetics paper describes
SS-H12, a different shared-sweep statistic, and is not the LASSI method paper.

### Completed LD-decay audit and primary MAF decision

The CT LD-decay audit completed under both prespecified frequency thresholds.
With MAF>=0.05, 4,033,512 valid anchor-neighbor LD pairs were summarized.
Mean r2 in the 0.5-1.5-kb baseline interval was 0.41449, giving a one-third
threshold of 0.13816. The first 1-kb bin below this threshold was centered at
55.5 kb; the five-consecutive-bin stability diagnostic began at 57.5 kb.

The MAF>=0.01 sensitivity included 7,272,074 valid pairs. Its baseline mean r2
was 0.18648 and its direct/stable crossings were 63.5/67.5 kb. The sensitivity
therefore changes the absolute estimate modestly but not its order of magnitude.

The production decision is to calibrate the LASSI window from the MAF>=0.05
LD curve. This is because the planned LASSI scan uses the same MAF>=0.05 marker
panel, so the SNP-delimited window should be calibrated from that panel rather
than from a broader marker set. MAF>=0.01 remains a documented sensitivity and
is not averaged with the primary estimate. This choice is tied to the LASSI
marker panel rather than to an intrinsic phasing requirement.

The selected physical decay estimate is therefore 55.5 kb according to the
literature-matched first-crossing rule, with 57.5 kb retained as the stability
diagnostic. The final LASSI SNP-window count remains to be obtained by
translating 55.5 kb to the empirical marker count in the frozen MAF>=0.05
LASSI panel after visual review of the static and interactive curves.

### Cross-analysis figure policy

Scientific figures are now accompanied, where useful, by self-contained HTML
versions for interactive inspection. Individual-level hover text reports sample
ID, population and plotted quantities; pairwise KING figures report both sample
IDs; LD-decay HTML reports distance bin, mean r2, pair count and uncertainty.
Static PNG/PDF files remain the manuscript-oriented outputs. Captions are
versioned in config/figure_registry.tsv and docs/figure_captions.md.

### LASSI window units: SNP count, not physical distance

The 55.5-kb LD-decay estimate is an intermediate calibration scale only.
LASSI itself uses an SNP-delimited sliding window, so the production parameter
must be expressed as a number of SNPs rather than kb.

The workflow therefore converts the selected 55.5-kb physical interval to an
empirical SNP count using the same final CT MAF>=0.05 marker panel used for
the LD audit and planned LASSI scan. At each of the 2,200 deterministic
genome-wide LD anchors, the number of retained SNPs in a centered window of
total width 55.5 kb is counted. The rounded genome-wide median of those counts
is the candidate LASSI winsize. The 57.5-kb stable LD crossing is converted in
parallel as a robustness diagnostic. Winstep is then set to approximately 10%
of winsize, rounded to an integer.

This avoids substituting physical distance directly into LASSI and avoids
using a single genome-wide average marker density in the presence of local
variation in WGS SNP density. The final winsize and winstep remain unfrozen
until the empirical count distribution has been inspected.


### Final LASSI comparison populations and population-specific windows

The comparative LASSI design is restricted to four populations with distinct
roles. Cinque Terre (CT) is the focal discovery population. TSI and IBS are
Southern-European comparator populations chosen for geographic and genetic
proximity. CEU is retained as an empirical LASSI benchmark / validation
reference because Harris and DeGiorgio applied LASSI to 1000 Genomes CEU and
reported established European sweep candidates, including the LCT/MCM6 region.

The original LASSI human application used a common 117-SNP window for CEU and
YRI. The present study deliberately does not copy that choice across all
populations. Instead, CT, CEU, TSI and IBS each receive a population-specific
MAF>=0.05 LD-decay analysis and an independently derived SNP-delimited LASSI
window. This is an explicit methodological extension: each population's
analysis window is calibrated to its own background LD scale.

For CT, calibration is complete. Across 2,200 deterministic anchors, 55.5-kb
windows contain a median of 99 SNPs (IQR 68-136), while the 57.5-kb stability
interval contains a median of 103 SNPs (IQR 70-141). The small difference
supports winsize=99 SNP and winstep=10 SNP for CT.

The population-specific calibration is now complete. CEU crosses the one-third
baseline threshold at 53.5 kb and yields a median 116 SNPs per calibrated window
(IQR 78.75-156), giving winsize=116 and winstep=12. TSI crosses at 52.5 kb,
with median 112 SNPs (IQR 77-154), giving winsize=112 and winstep=11. IBS
crosses at 53.5 kb, with median 116 SNPs (IQR 78.75-156), giving winsize=116
and winstep=12. For all three 1000 Genomes populations, the direct first
crossing is already stable for at least five consecutive bins.

Across CT, CEU, TSI and IBS the selected physical scales are therefore tightly
clustered at 52.5-55.5 kb, whereas the SNP-delimited windows range from 99 to
116 SNPs. This distinction is biologically and technically useful: background
LD scale is similar at this resolution, while different population-specific
MAF>=0.05 marker densities translate that scale into different LASSI SNP
windows. The final scan parameters are CT 99/10, CEU 116/12, TSI 112/11 and
IBS 116/12 (winsize/winstep, in SNP).

The CEU result is also an empirical implementation benchmark: the independently
derived CEU winsize of 116 SNPs is extremely close to the 117-SNP window used
in the original human LASSI application. This agreement is treated as a
plausibility check, not as an exact replication target, because preprocessing,
filtering and calibration details differ.

The three 1000 Genomes comparison populations are already phased in the Phase 3
source panel, so they do not require re-phasing for LASSI. The phasing branch is
required for CT only. CEU benchmarking will be interpreted qualitatively and
positionally against the published LASSI CEU scan rather than as an expectation
of numerically identical T statistics.


### Full 1000 Genomes EUR LD-decay context

LD decay is additionally estimated in FIN and GBR with the same MAF>=0.05,
100-anchor-per-autosome and one-third-baseline protocol used for CEU, TSI and
IBS. This extension is descriptive rather than a change to the LASSI scan
design. It allows the CT curve to be interpreted against the complete set of
five 1000 Genomes EUR populations and avoids presenting only the populations
selected for downstream LASSI.

No LASSI window is derived or used for FIN or GBR in the primary analysis.
The frozen selection scans remain CT, CEU, TSI and IBS with their existing
population-specific windows. If FIN or GBR are later promoted to formal LASSI
comparators, their SNP-window conversion must then be performed explicitly
from their own LD decay rather than borrowing a window from another population.


### Completed full European LD-decay context

The full MAF>=0.05 LD-decay context is now complete for CT and all five
1000 Genomes EUR populations. FIN has baseline mean r2=0.41292, first crossing
55.5 kb and five-bin stable crossing 58.5 kb. GBR has baseline mean r2=0.42346,
with both first and stable crossing at 55.5 kb.

Together with CEU (53.5 kb), TSI (52.5 kb), IBS (53.5 kb) and CT (55.5 kb),
the first-crossing estimates span only 52.5-55.5 kb across all six populations.
This supports the conclusion that the physical LD scale is broadly similar
across the European context at the resolution used here. FIN is the only
population in which the five-bin stability diagnostic extends appreciably
beyond the first crossing (58.5 kb versus 55.5 kb), but this does not alter the
descriptive first-crossing comparison.

FIN and GBR remain contextual populations only for the selection branch.
Their LD curves are included in the static and interactive six-population
comparison, but no production LASSI winsize is derived because they are not
part of the planned LASSI scan. CT/CEU/TSI/IBS winsize/winstep parameters
remain frozen and unchanged.

### CT phasing: SHAPEIT2 reference-assisted production design

The CT production phasing method is SHAPEIT2 v2.r904 with the phased
1000 Genomes Phase 3 EUR reference panel. This is a deliberate small-study
reference-assisted design: the SHAPEIT2 documentation states that external
reference haplotypes are particularly useful when phasing fewer than about
100 study individuals, whereas CT contains 46.

The effective population size is fixed at Ne=11,418, the European/CEU value
recommended by the SHAPEIT2 documentation. The common-variant conditioning
window is fixed at 0.5 Mb because the documentation specifically advises this
value for sequence-derived genotypes rather than the 2-Mb GWAS default.
Other model/MCMC settings remain at the documented SHAPEIT2 defaults:
100 conditioning states, 7 burn-in, 8 pruning and 20 main iterations.
The --no-mcmc shortcut is not used because it is recommended only when the
study contains typically fewer than 10 individuals.

Reference alignment is checked formally with SHAPEIT2 -check. In
reference-assisted mode, a study SNP must also occur in the reference with
compatible alleles; study-only or incompatible positions are reported in the
alignment exclusion file and are removed from the phasing command with
--exclude-snp. The chromosome-wise bcftools overlap audit remains as an
independent quantitative summary of marker retention, but the actual SHAPEIT2
-check output is authoritative for the production exclusion list.

The reference scope is not an open parameter: it is 1000 Genomes EUR. The local HAP/LEGEND/SAMPLE files are the complete 2503-sample Phase 3 reference, and the EUR restriction is implemented with SHAPEIT2 `--include-grp` using a one-line group file containing `EUR`. This is preferable to duplicating/subsetting the large reference files on disk and is directly supported by SHAPEIT2. GRCh37
genetic maps are supplied explicitly. A fixed seed and one SHAPEIT2 thread per
chromosome are used for exact reproducibility, while Snakemake may parallelize
different chromosomes.


### Completed SHAPEIT2 EUR overlap preflight

The chromosome-wise CT-versus-1000 Genomes EUR overlap preflight completed
successfully across all 22 autosomes. The CT MAF>=0.05 panel reconstructed
directly from the final QC VCF contains 5,007,326 SNPs, exactly matching the
existing PLINK panel on every chromosome.

After subsetting the phased 1000 Genomes Phase 3 data to the 503 EUR samples
and retaining variants polymorphic within EUR, 4,883,824 CT SNPs have an exact
CHR:POS:REF:ALT match. This is 97.5336% of the CT common-variant panel.
There are 4,884,025 same-position overlaps, of which only 201 are
allele-discordant. Thus 123,502 CT MAF>=0.05 SNPs lack an exact match in the
polymorphic-EUR preflight panel.

Coverage is consistently high across chromosomes: the lowest exact-match
fraction is 95.1719% on chromosome 21 and the highest is 98.1391% on
chromosome 2. The result supports use of the fixed 1000 Genomes EUR reference
for CT phasing.

This overlap audit is intentionally not treated as the final exclusion list.
It counts only variants that remain polymorphic after EUR subsetting, whereas
the actual SHAPEIT2 HAP/LEGEND/SAMPLE reference representation may differ in
site content. The authoritative production exclusion list is therefore the
one generated by SHAPEIT2 `-check`.

Operational correction: standard genotype phasing with an external reference
uses SHAPEIT2's default `-phase` mode. The `-assemble` mode is specifically
for phase-informative reads (PIRs) and is not used in the present analysis.
The production sequence is therefore: prepare CT MAF>=0.05 chromosome VCF ->
SHAPEIT2 `-check` with EUR reference and GRCh37 map -> inspect the generated
`.snp.strand.exclude` list -> standard reference-assisted SHAPEIT2 phasing
with that exclusion list.


### Historical SHAPEIT2 run versus current EUR-only reference design

Inspection of the archived 2025 CT phasing logs shows that the previous run
used the complete `1000GP_Phase3_chr*.hap.gz`, `.legend.gz`, and
`1000GP_Phase3.sample` resources. The SAMPLE file contains 2503 individuals,
so that historical run was not physically restricted to EUR and its logged
command lines do not show an explicit `--include-grp EUR` filter.

The current production pipeline intentionally differs at this point: it uses
the same complete Phase 3 files but supplies `--include-grp` with the group
identifier `EUR`, yielding the intended 503-sample European reference
without generating duplicated HAP/LEGEND/SAMPLE resources. The same group
filter is applied during both SHAPEIT2 `-check` and production phasing so that
alignment diagnostics and phasing use the identical reference subset.


### SHAPEIT2 two-stage alignment check

SHAPEIT2 reference-assisted phasing requires study variants to occur in the
reference with compatible alleles. Its `-check` mode reports study-only and
allele-incompatible sites in `.snp.strand` and writes their physical
positions to `.snp.strand.exclude`. When such problems are found, SHAPEIT2
terminates the diagnostic command non-zero; this is an expected control-flow
condition rather than evidence that the reference panel itself is invalid.

The workflow therefore implements the documented two-stage procedure. First,
an unfiltered `-check` is run against the 1000G EUR subset. A non-zero return
is accepted only if both diagnostic files are non-empty. Second, the generated
exclusion list is passed back to SHAPEIT2 with `--exclude-snp` in a new
`-check`. This second check must return zero. Production phasing has an
explicit dependency on that per-chromosome PASS sentinel.

The first observed chromosome under this workflow was chromosome 12:
240,237 CT MAF>=0.05 SNPs were read, the EUR group filter correctly retained
1,006 reference haplotypes (503 diploid individuals) and excluded 4,002
non-EUR haplotypes, and SHAPEIT2 identified 3,333 study SNPs missing from the
EUR reference plus 143 misaligned sites. These 3,476 positions are exactly the
kind of variants the formal exclusion step is intended to remove before the
clean second check and subsequent phasing.


### SHAPEIT2 conditioning-window choice (`--window 0.5`)

The production CT phasing uses a 0.5-Mb SHAPEIT2 conditioning window. This is
not derived from the LASSI LD-decay analysis and should not be interpreted as
an estimate of the physical extent of LD in Cinque Terre. It is a model
parameter controlling the local genomic interval over which SHAPEIT2 selects
conditioning haplotypes.

Primary software documentation:
SHAPEIT2 documentation, "Model parameters: Window size W (--window)".
The documentation states that the default is approximately 2 Mb for GWAS
datasets, but that the developers' experiments suggest 0.5 Mb may give better
results for sequence data; the option table likewise advises 0.5 Mb for
genotypes derived from sequencing.

Published methodological support:

Delaneau O, Howie B, Cox AJ, Zagury J-F, Marchini J. 2013.
**Haplotype Estimation Using Sequencing Reads.**
*American Journal of Human Genetics* 93(4):687-696.
doi:10.1016/j.ajhg.2013.09.002. PMID:24094745. PMCID:PMC3791270.

Role in this project:
- directly studies SHAPEIT2 in sequencing-derived genotype data;
- includes high-coverage Illumina samples together with European 1000 Genomes
  reference haplotypes/genotypes;
- in the standard SHAPEIT2 experiment without phase-informative reads, uses
  W=0.5 Mb and K=100;
- therefore provides a close methodological precedent for a 0.5-Mb
  conditioning window in sequence-derived human genotypes with European
  reference support.

Sharp K, Kretzschmar W, Delaneau O, Marchini J. 2016.
**Phasing for medical sequencing using rare variants and large haplotype
reference panels.**
*Bioinformatics* 32(13):1974-1980.
doi:10.1093/bioinformatics/btw065. PMID:27153703. PMCID:PMC4920110.

Role in this project:
- explicitly states that SHAPEIT2 had previously shown good performance with
  a 0.5-Mb window for unphased genotypes derived from sequencing, citing the
  earlier SHAPEIT2 work and its Supplementary Figure S3;
- consequently uses 0.5 Mb throughout its reference-based sequencing phasing
  experiments;
- provides an independent later methodological use of the same window in the
  context of sequencing-derived genotypes and external haplotype reference
  panels.

Project interpretation:
- 0.5 Mb is adopted because it is the developers' sequence-data recommendation
  and is supported by published SHAPEIT2 sequencing/reference-panel studies.
- The choice is not tuned to the observed CT LD curve and is therefore
  independent of the approximately 55-kb physical LD scale used only to
  calibrate the SNP-count window for LASSI.
- We do not claim that 0.5 Mb is a universally optimal phasing window for all
  WGS datasets; rather, it is a literature-supported, software-recommended
  setting that matches the data type and reference-assisted design used here.


### Completed CT SHAPEIT2 alignment audit

The formal two-stage SHAPEIT2 alignment audit has completed on all 22
autosomes using the 1000 Genomes Phase 3 reference restricted internally to
EUR (503 individuals; 1006 reference haplotypes).

Genome-wide results:
- CT MAF>=0.05 SNPs entering the audit: 5,007,326;
- missing from the EUR reference: 90,392;
- allele-misaligned between CT and EUR reference: 2,651;
- actual SNPs excluded by SHAPEIT2 before phasing: 93,043;
- SNPs retained for production phasing: 4,914,283;
- retained fraction: 0.9814186254 (98.14%);
- minimum chromosome-specific retained fraction: 0.9578043165 (chr21);
- maximum chromosome-specific retained fraction: 0.9873842488 (chr4);
- all 22 post-exclusion SHAPEIT2 `-check` runs completed successfully.

The `.snp.strand.exclude` files contain 94,923 non-empty lines in total,
which is slightly larger than the 93,043 variants actually excluded. The
version-controlled audit therefore records SHAPEIT2's reported numbers of
included/excluded study SNPs from the second check logs, rather than using raw
line counts from the exclusion files as the authoritative variant count.

This audit is considered passed and authorizes parameter benchmarking and
subsequent production phasing.


### SHAPEIT2 chr20 K=400 runtime benchmark

Before changing the production conditioning-state count, the workflow performs
a dedicated chromosome-20 benchmark with `--states 400`, `--window 0.5`,
`--effective-size 11418`, the standard 7/8/20 burn/prune/main schedule, a
fixed seed, and one SHAPEIT2 thread. Chromosome 20 is used because the formal
alignment audit leaves 111,416 CT SNPs for phasing there, providing a
representative but tractable chromosome for measuring computational cost.

The benchmark is deliberately isolated from production outputs. Its purpose is
to measure wall time and memory under K=400 and to confirm that the proposed
accuracy-oriented setting is computationally practical. It is not used to
select a parameter based on downstream LASSI results. Production remains at
the previously configured state count until this benchmark is reviewed.

The benchmark uses Snakemake's native benchmark recording and writes its
phasing output under `results/selection/phasing/shapeit2/benchmark/`, with
resource metrics under `benchmarks/selection/phasing/`.

Observed result: the chr20 K=400, W=0.5 Mb, Ne=11,418, 7/8/20 MCMC,
thread=1 run completed successfully on 111,416 retained CT SNPs in exactly
2,400 seconds (40.0 minutes) according to SHAPEIT2. The log confirmed
`400 states per window [400 H + 0 PM + 0 R + 0 COV]`, 35 MCMC iterations,
and the fixed seed 15052011. This demonstrates that K=400 is computationally
practical at chromosome scale on the current workstation. The Snakemake benchmark recorded max RSS 836.42 MB, max USS 832.88 MB,
max PSS 833.02 MB, mean load 99.49%, and CPU time 2,387.10 s for 2,399.23 s
wall time. Thus the single-thread run is strongly CPU-bound and uses less than
1 GB resident memory. This supports keeping one SHAPEIT2 thread per chromosome
for deterministic reproducibility while exploiting machine-level parallelism
by running independent chromosomes concurrently through Snakemake. Practical
concurrency should be chosen from the workstation's available CPU cores,
memory, and thermal behavior rather than by increasing SHAPEIT2 threads within
a chromosome.


### Final CT SHAPEIT2 conditioning-state decision

Following the completed chromosome-20 benchmark, the production CT phasing
state count is frozen at `--states 400`. The benchmark phased 111,416 chr20
SNPs in 2,400 s with one thread, max RSS 836.42 MB, and mean CPU load 99.49%.
This demonstrates that K=400 is computationally practical on the available
workstation while providing a substantially larger conditioning set than the
SHAPEIT2 default K=100. The production configuration therefore uses 400
Hamming-distance-selected conditioning haplotypes per window, with no
additional perfect-match, random, or coverage-based states.

Genome-wide execution remains one SHAPEIT2 thread per chromosome for exact
seed reproducibility. Parallelism is delegated to Snakemake by running
independent chromosomes concurrently, avoiding within-chromosome
multithreading while efficiently using the workstation's available CPU cores.


### Completed CT SHAPEIT2 production phasing

Genome-wide CT phasing completed successfully across all 22 autosomes using
SHAPEIT2 v2.r904 with the frozen production settings: 1000 Genomes Phase 3 EUR
reference (1,006 reference haplotypes), Ne=11,418, W=0.5 Mb, K=400 Hamming-
selected conditioning states, 7 burn-in, 8 pruning and 20 main MCMC iterations,
one SHAPEIT2 thread per chromosome, and seed 15052011.

The 22 chromosome jobs were executed concurrently through Snakemake with a
global `--cores 10` limit. The aggregate target completed 23/23 workflow
steps successfully on 2026-09-24. Total wall-clock time for the genome-wide
production invocation was 12:08:31.076032. All expected chromosome-level
HAP/SAMPLE outputs and SHAPEIT2 phase logs were present at completion.

This completion establishes the phased CT common-variant panel as the source
for downstream LASSI input construction. No additional re-phasing is required
before LASSI preprocessing unless a later sensitivity analysis is explicitly
introduced.


## Post-phasing SHAPEIT2 production audit

Before the phased CT panel is accepted as the production input for LASSI, the
workflow performs an exhaustive internal-consistency audit across all 22
autosomes. The audit verifies gzip/CRC integrity, 46-sample identity and order
against each chromosome-specific input VCF, 92 binary haplotype states per
site, chromosome/position ordering, absence of duplicate phased sites and
physical positions, and exact per-chromosome and genome-wide site counts
against the completed two-stage SHAPEIT2 alignment audit.

Genotype preservation is checked exhaustively rather than by spot sampling:
for every phased site, the two SHAPEIT2 haplotypes are collapsed back to an
unphased dosage and compared with every non-missing source VCF genotype.
Allele swaps and non-palindromic strand complements are handled explicitly.
Any dosage mismatch is treated as an audit failure.

The chromosome-20 K=400 benchmark and chromosome-20 production run used the
same study input, EUR reference restriction, map, seed, one-thread execution,
window, effective population size and MCMC settings. Their decompressed HAPS
content and SAMPLE files are therefore compared as an empirical repeatability
test of the frozen single-thread seeded workflow.

The audit also checks each production phase log for completed main iterations,
graph normalization, haplotype solving, the expected SNP/sample/reference
counts and the K=400 model, while separately recording whether seed, thread,
0.5-Mb window and Ne=11,418 are parsed from every log.

These checks establish technical integrity, genotype conservation and
repeatability. They do not estimate switch error or statistical phasing
accuracy, which would require independent phase truth (for example trios or a
truth set) or a prespecified external-method sensitivity analysis.


### Completed post-phasing internal-consistency audit

The production CT SHAPEIT2 panel passed the post-phasing audit across all 22
autosomes. The observed total of 4,914,283 phased SNPs exactly matched the
count retained by the formal SHAPEIT2 alignment audit. All chromosomes
contained 46 study samples and 92 haplotypes per site. No site-identity
mismatches, non-binary haplotype values, duplicate exact sites, duplicate
physical positions, or out-of-order coordinates were detected.

Genotype conservation was checked exhaustively at all phased sites. A total of
226,032,551 non-missing source VCF genotype calls were compared with diploid
dosages reconstructed from the two SHAPEIT2 haplotypes; zero dosage mismatches
were observed. A further 24,467 source genotype calls were missing and were
therefore not used for genotype-conservation comparison.

All core production-log checks passed, including parsing of seed 15052011,
one thread per chromosome, the 0.5-Mb window and Ne=11,418.

Chromosome 20 provided an empirical repeatability test. The production and
states-400 benchmark HAPS files had identical decompressed SHA256 content and
the SAMPLE files were also identical, giving a repeatability PASS for the
fixed-seed, single-thread production configuration.

The phased CT panel is structurally ready for LASSI input construction. These
checks validate integrity, genotype conservation and repeatability; they do
not estimate switch error or absolute phase accuracy in the absence of
independent phase truth.
