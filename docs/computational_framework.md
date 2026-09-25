# Computational framework

This document summarizes the main software, utilities, and computational tools used throughout the pipeline, together with their primary purposes and reference citations.

---

# Analytical tools

## PLINK

Purpose:
Population genomics analyses, variant-level filtering, LD pruning, IBS estimation, ROH detection, and relatedness exploration.

Main functionalities used:

* variant filtering
* LD pruning
* IBS / pairwise similarity estimation
* PLINK `--genome` checks
* ROH analyses

Citations:

Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A. R., Bender, D., Maller, J., Sklar, P., de Bakker, P. I. W., Daly, M. J., & Sham, P. C. (2007). PLINK: A tool set for whole-genome association and population-based linkage analyses. *American Journal of Human Genetics, 81*(3), 559–575. https://doi.org/10.1086/519795

Chang, C. C., Chow, C. C., Tellier, L. C. A. M., Vattikuti, S., Purcell, S. M., & Lee, J. J. (2015). Second-generation PLINK: Rising to the challenge of larger and richer datasets. *GigaScience, 4*(1), 7. https://doi.org/10.1186/s13742-015-0047-8

---

## KING

Purpose:
Robust relatedness inference and kinship estimation.

Main functionalities used:

* relatedness estimation
* duplicate detection
* close-relative checks

In this workflow, KING is mainly used to support relatedness exploration and to check for unexpected duplicates or close relatives.

Citation:

Manichaikul, A., Mychaleckyj, J. C., Rich, S. S., Daly, K., Sale, M., & Chen, W. M. (2010). Robust relationship inference in genome-wide association studies. *Bioinformatics, 26*(22), 2867–2873. https://doi.org/10.1093/bioinformatics/btq559

---

## EIGENSOFT

Purpose:
Population structure inference and principal component analysis.

Main functionalities used:

* PCA
* smartpca

Citation:

Patterson, N., Price, A. L., & Reich, D. (2006). Population structure and eigenanalysis. *PLOS Genetics, 2*(12), e190. https://doi.org/10.1371/journal.pgen.0020190

---

## SHAPEIT2

Purpose:
Haplotype phasing and inference of chromosomal phase.

Main functionalities used:

* reference-supported phasing
* haplotype reconstruction
* preparation of phased data for LASSI and optional IBD-based analyses

Citation:

Delaneau, O., Zagury, J. F., & Marchini, J. (2013). Improved whole-chromosome phasing for disease and population genetic studies. *Nature Methods, 10*(1), 5–6. https://doi.org/10.1038/nmeth.2307

---

## lassip: saltiLASSI and original LASSI

Purpose:
Likelihood-based haplotype-frequency-spectrum scans for selective-sweep candidates.

Pinned implementation:

* `lassip` v1.2.1
* commit `a6a9d18c2323330fbf74d5a490f9e9c4ebe41d7c`

Primary functionality used:

* saltiLASSI spatial composite-likelihood statistic Lambda (`--salti`)

Secondary functionality used:

* original LASSI likelihood-ratio statistic T (`--lassi`)

The production model uses phased input, K=10 and Harris & DeGiorgio Model D; in lassip v1.2.1 this maps to `--lassi-choice 4`.

The scan identifies candidate/outlier regions for downstream interpretation; it does not by itself prove adaptation or identify causal variants.

Citations:

Harris, A. M., & DeGiorgio, M. (2020). A likelihood approach for uncovering selective sweep signatures from haplotype data. *Molecular Biology and Evolution, 37*(10), 3023-3046. https://doi.org/10.1093/molbev/msaa115

DeGiorgio, M., & Szpiech, Z. A. (2022). A spatially aware likelihood test to detect sweeps from haplotype distributions. *PLOS Genetics, 18*, e1010134. https://doi.org/10.1371/journal.pgen.1010134

---

## hap-IBD

Purpose:
Detection of identity-by-descent haplotype segments.

Main functionalities used:

* IBD segment detection
* optional preparation of IBD segment inputs for recent demographic inference

This tool is considered part of the optional demographic branch.

Citation:

Zhou, Y., Browning, S. R., & Browning, B. L. (2020). A fast and simple method for detecting identity-by-descent segments in large-scale data. *American Journal of Human Genetics, 106*(4), 426–437. https://doi.org/10.1016/j.ajhg.2020.02.010

---

## IBDNe

Purpose:
Inference of recent effective population size from IBD segment distributions.

Main functionalities used:

* recent demographic inference
* effective population size estimation
* Ne trajectory reconstruction

This tool is considered part of the optional demographic branch.

Citations:

Browning, S. R., & Browning, B. L. (2015). Accurate non-parametric estimation of recent effective population size from segments of identity by descent. *American Journal of Human Genetics, 97*(3), 404–418. https://doi.org/10.1016/j.ajhg.2015.07.012

Fenner, J. N. (2005). Cross-cultural estimation of the human generation interval for use in genetics-based population divergence studies. *American Journal of Physical Anthropology, 128*(2), 415–423. https://doi.org/10.1002/ajpa.20188

---

## g:Profiler

Purpose:
Exploratory functional enrichment and over-representation analysis.

Main functionalities used:

* gene ontology enrichment
* pathway enrichment
* functional annotation of genes overlapping or near candidate regions

In this workflow, g:Profiler is used only for exploratory biological contextualization, not as independent validation of selection.

Citation:

Raudvere, U., Kolberg, L., Kuzmin, I., Arak, T., Adler, P., Peterson, H., & Vilo, J. (2019). g:Profiler: A web server for functional enrichment analysis and conversions of gene lists (2019 update). *Nucleic Acids Research, 47*(W1), W191–W198. https://doi.org/10.1093/nar/gkz369

---

## STRING

Purpose:
Protein-protein interaction network analysis and functional interaction exploration.

Main functionalities used:

* functional interaction networks
* protein interaction visualization
* exploratory biological contextualization of candidate-region genes

In this workflow, STRING is used only as a downstream exploratory interpretation tool.

Citation:

Szklarczyk, D., Gable, A. L., Nastou, K. C., Lyon, D., Kirsch, R., Pyysalo, S., Doncheva, N. T., Legeay, M., Fang, T., Bork, P., Jensen, L. J., & von Mering, C. (2021). The STRING database in 2021: Customizable protein-protein networks, and functional characterization of user-uploaded gene/measurement sets. *Nucleic Acids Research, 49*(D1), D605–D612. https://doi.org/10.1093/nar/gkaa1074

---

# Workflow and infrastructure

## Snakemake

Purpose:
Workflow management and reproducible pipeline orchestration.

Main functionalities used:

* rule-based workflow execution
* dependency management
* reproducible analyses
* modular pipeline development

Citation:

Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., & Köster, J. (2021). Sustainable data analysis with Snakemake. *F1000Research, 10*, 33. https://doi.org/10.12688/f1000research.29032.2

---

## bcftools

Purpose:
Manipulation, normalization, filtering, and harmonization of VCF/BCF files.

Main functionalities used:

* variant filtering
* dataset harmonization
* VCF manipulation
* dataset merging
* variant normalization

Citation:

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience, 10*(2), giab008. https://doi.org/10.1093/gigascience/giab008

---

## samtools

Purpose:
Manipulation and indexing of sequencing and alignment-related files.

Main functionalities used:

* FASTA indexing
* alignment inspection, if needed
* sequence handling

Citation:

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience, 10*(2), giab008. https://doi.org/10.1093/gigascience/giab008

---

## tabix

Purpose:
Indexing and rapid genomic querying of compressed genomic files.

Main functionalities used:

* indexing compressed VCF files
* genomic interval queries

Citation:

Li, H. (2011). Tabix: Fast retrieval of sequence features from generic TAB-delimited files. *Bioinformatics, 27*(5), 718–719. https://doi.org/10.1093/bioinformatics/btq671

---

## Conda

Purpose:
Environment and dependency management.

Main functionalities used:

* reproducible software environments
* package management
* workflow-specific dependency isolation

Citation:

Anaconda Software Distribution. (2020). *Anaconda documentation*. Anaconda Inc. https://docs.anaconda.com/

---

# Utility and shell tools

## GNU awk

Purpose:
Text processing, parsing, and manipulation of tabular outputs.

Main functionalities used:

* tabular output parsing
* field extraction
* lightweight text-based transformations

Citation:

Aho, A. V., Kernighan, B. W., & Weinberger, P. J. (2023). *The AWK programming language* (2nd ed.). Addison-Wesley.

---

## GNU coreutils

Purpose:
General shell-level file and text manipulation utilities.

Main utilities potentially used:

* cat
* cp
* mv
* rm
* sort
* uniq
* wc
* cut
* paste
* tr

Citation:

Free Software Foundation. *GNU Core Utilities*. https://www.gnu.org/software/coreutils/

---

## GNU grep

Purpose:
Pattern matching and text searching.

Main functionalities used:

* filtering outputs
* searching log files
* pattern extraction

Citation:

Free Software Foundation. *GNU Grep manual*. https://www.gnu.org/software/grep/

---

## GNU sed

Purpose:
Stream editing and text transformation.

Main functionalities used:

* text substitution
* pipeline text manipulation

Citation:

Free Software Foundation. *GNU sed manual*. https://www.gnu.org/software/sed/

---

## Bash

Purpose:
Shell scripting and command-line workflow execution.

Main functionalities used:

* pipeline scripting
* workflow orchestration
* command chaining

Citation:

Free Software Foundation. *Bash reference manual*. https://www.gnu.org/software/bash/
