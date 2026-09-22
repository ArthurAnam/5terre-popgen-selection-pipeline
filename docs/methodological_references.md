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

#### Frozen production ROH definition and FROH

The production call uses the joint MAF>=0.05, unpruned marker panel with a
50-SNP sliding window, a minimum of 50 SNP per called ROH, a minimum physical
length of 1.5 Mb, density <=50 kb/SNP, gap <=500 kb, <=5 missing calls per
window, <=1 heterozygous call per window, and window hit threshold 0.05.

The distinction between the two SNP-count parameters is explicit:
`--homozyg-window-snp 50` defines the size of the moving scan window, whereas
`--homozyg-snp 50` requires each final called ROH to contain at least 50 SNPs.
In this dense WGS-derived panel, the 1.5-Mb physical threshold is the dominant
biological scale and the 50-SNP final-run minimum is mainly a technical floor.

PLINK distinguishes the local scanning-window heterozygote limit
(`--homozyg-window-het`) from the whole-run heterozygote cap
(`--homozyg-het`). The production analysis leaves `--homozyg-het` unset:
a fixed whole-run cap would become increasingly stringent as segment length
increases, while one heterozygote per 50-SNP scanning window provides local
tolerance for isolated genotype errors.

Primary FROH is the summed physical length of all autosomal ROH >=1.5 Mb divided
by 2.77e9 bp. This 2.77-Gb denominator is the total SNP-mappable autosomal
distance used in prior human FROH studies. The project's ~2.794-Gb
terminal-marker span is not used because it was defined for marker-density
auditing and bridges marker-free intervals rather than estimating a callable
autosomal denominator.

Production summaries retain continuous ROH lengths and additionally report
1.5-<5 Mb and >=5 Mb classes. These bins are descriptive: the primary phenotype
remains FROH based on all ROH >=1.5 Mb.

#### Completed production ROH and FROH

The frozen production workflow completed on all 549 individuals and yielded
2,885 ROH >=1.5 Mb, of which 234 were >=5 Mb. The segment count is internally
consistent with the summed population-level N_ROH values.

CT has the highest mean primary FROH (0.01321; mean total ROH 36.58 Mb), followed
by FIN (0.00701; 19.41 Mb), IBS (0.00492; 13.63 Mb), GBR (0.00431; 11.94 Mb),
CEU (0.00311; 8.62 Mb) and TSI (0.00294; 8.14 Mb). However, CT median FROH
(0.00734) is very close to FIN (0.00703), so the CT-vs-FIN contrast is not
well described as a simple uniform upward shift.

The length composition is more distinctive. CT has a median of one ROH >=5 Mb
and mean long-ROH burden of 18.36 Mb, whereas every reference population has a
median of zero >=5-Mb ROH. On average, ROH >=5 Mb account for ~50.2% of the CT
total >=1.5-Mb burden, compared with ~12.0% in FIN, 14.1% in CEU, 20.3% in GBR,
22.9% in TSI and 27.6% in IBS. This pattern is consistent with a stronger
recent-autozygosity component in at least a subset of CT individuals, but that
interpretation must be based on the individual distribution rather than the
population mean alone.

The next analysis therefore treats FROH>=1.5 Mb as the primary endpoint and
FROH>=5 Mb as a secondary length-specific endpoint. CT is compared separately
with each EUR reference population using two-sided Mann-Whitney tests, Holm
correction across the five prespecified CT-vs-reference contrasts, and Cliff's
delta as a distributional effect-size measure. Means, medians, quartiles and
individual scatter are reported alongside the tests.


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
