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
