# Computational framework

This document summarizes the main software, utilities, and computational tools used throughout the pipeline, together with their primary purposes and reference citations.

---

# Analytical tools

## PLINK

Purpose:
Population genomics analyses, variant-level QC, LD pruning, IBS estimation, ROH detection, and relatedness filtering.

Main functionalities used:
- QC filtering
- LD pruning
- IBS/IBD estimation
- ROH analyses

Citations:

Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A. R., Bender, D., Maller, J., Sklar, P., de Bakker, P. I. W., Daly, M. J., & Sham, P. C. (2007). PLINK: A tool set for whole-genome association and population-based linkage analyses. American Journal of Human Genetics, 81(3), 559–575. https://doi.org/10.1086/519795

Chang, C. C., Chow, C. C., Tellier, L. C. A. M., Vattikuti, S., Purcell, S. M., & Lee, J. J. (2015). Second-generation PLINK: Rising to the challenge of larger and richer datasets. GigaScience, 4(1), 7. https://doi.org/10.1186/s13742-015-0047-8

---

## KING

Purpose:
Robust relatedness inference and kinship estimation.

Main functionalities used:
- Relatedness detection
- Duplicate identification

Citation:

Manichaikul, A., Mychaleckyj, J. C., Rich, S. S., Daly, K., Sale, M., & Chen, W. M. (2010). Robust relationship inference in genome-wide association studies. Bioinformatics, 26(22), 2867–2873. https://doi.org/10.1093/bioinformatics/btq559

---

## EIGENSOFT

Purpose:
Population structure inference and principal component analysis.

Main functionalities used:
- PCA
- smartpca

Citation:

Patterson, N., Price, A. L., & Reich, D. (2006). Population structure and eigenanalysis. PLOS Genetics, 2(12), e190. https://doi.org/10.1371/journal.pgen.0020190

---

## SHAPEIT2

Purpose:
Haplotype phasing and inference of chromosomal phase.

Main functionalities used:
- Reference-guided phasing
- Haplotype reconstruction

Citation:

Delaneau, O., Zagury, J. F., & Marchini, J. (2013). Improved whole-chromosome phasing for disease and population genetic studies. Nature Methods, 10(1), 5–6. https://doi.org/10.1038/nmeth.2307

---

## LASSI

Purpose:
Likelihood-based detection of selective sweep regions from haplotype data.

Main functionalities used:
- Sweep likelihood estimation
- Hard/soft sweep inference (haplotype/s under selective sweep process)

Citation:

Harris, A. M., & DeGiorgio, M. (2020). Identifying and classifying shared selective sweeps from multilocus data. Genetics, 215(1), 143–171. https://doi.org/10.1534/genetics.120.303049

---

## hap-IBD

Purpose:
Detection of identity-by-descent haplotype segments.

Main functionalities used:
- IBD segment detection

Citation:

Zhou, Y., Browning, S. R., & Browning, B. L. (2020). A fast and simple method for detecting identity-by-descent segments in large-scale data. American Journal of Human Genetics, 106(4), 426–437. https://doi.org/10.1016/j.ajhg.2020.02.010

---

## IBDNe

Purpose:
Inference of recent effective population size from IBD segment distributions.

Main functionalities used:
- Recent demographic inference
- Effective population size estimation

Citations:

Browning, S. R., & Browning, B. L. (2015). Accurate non-parametric estimation of recent effective population size from segments of identity by descent. American Journal of Human Genetics, 97(3), 404–418. https://doi.org/10.1016/j.ajhg.2015.07.012

Fenner, J. N. (2005). Cross-cultural estimation of the human generation interval for use in genetics-based population divergence studies. American Journal of Physical Anthropology, 128(2), 415–423. https://doi.org/10.1002/ajpa.20188

---

## g:Profiler

Purpose:
Functional enrichment and over-representation analysis (ORA).

Main functionalities used:
- Gene ontology enrichment
- Pathway enrichment
- Functional annotation

Citation:

Raudvere, U., Kolberg, L., Kuzmin, I., Arak, T., Adler, P., Peterson, H., & Vilo, J. (2019). g:Profiler: A web server for functional enrichment analysis and conversions of gene lists (2019 update). Nucleic Acids Research, 47(W1), W191–W198. https://doi.org/10.1093/nar/gkz369

---

## STRING

Purpose:
Protein-protein interaction network analysis and functional interaction exploration.

Main functionalities used:
- Functional interaction networks
- Protein interaction visualization

Citation:

Szklarczyk, D., Gable, A. L., Nastou, K. C., Lyon, D., Kirsch, R., Pyysalo, S., Doncheva, N. T., Legeay, M., Fang, T., Bork, P., Jensen, L. J., & von Mering, C. (2021). The STRING database in 2021: Customizable protein–protein networks, and functional characterization of user-uploaded gene/measurement sets. Nucleic Acids Research, 49(D1), D605–D612. https://doi.org/10.1093/nar/gkaa1074

---

# Workflow & infrastructure

## Snakemake

Purpose:
Workflow management and reproducible pipeline orchestration.

Main functionalities used:
- Rule-based workflow execution
- Dependency management
- Reproducible analyses

Citation:

Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., & Köster, J. (2021). Sustainable data analysis with Snakemake. F1000Research, 10, 33. https://doi.org/10.12688/f1000research.29032.2

---

## bcftools

Purpose:
Manipulation, normalization, filtering, and harmonization of VCF/BCF files.

Main functionalities used:
- Variant filtering
- Dataset harmonization
- VCF manipulation
- Dataset merging

Citation:

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

---

## samtools

Purpose:
Manipulation and indexing of sequencing and alignment files.

Main functionalities used:
- FASTA indexing
- Alignment inspection
- Sequence handling

Citation:

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

---

## tabix

Purpose:
Indexing and rapid genomic querying of compressed genomic files.

Main functionalities used:
- Indexing compressed VCF files
- Genomic interval queries

Citation:

Li, H. (2011). Tabix: Fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics, 27(5), 718–719. https://doi.org/10.1093/bioinformatics/btq671

---

## Conda

Purpose:
Environment and dependency management.

Main functionalities used:
- Reproducible software environments
- Package management

Citation:

Anaconda Software Distribution. (2020). Anaconda Documentation. Anaconda Inc. https://docs.anaconda.com/

---

# Utility & shell tools

## GNU awk

Purpose:
Text processing, parsing, and manipulation of tabular outputs.

Citation:

Aho, A. V., Kernighan, B. W., & Weinberger, P. J. (2023). The AWK programming language (2nd ed.). Addison-Wesley.

---

## GNU coreutils

Purpose:
General shell-level file and text manipulation utilities.

Main utilities potentially used:
- cat
- cp
- mv
- rm
- sort
- uniq
- wc
- cut
- paste
- tr

Citation:

Free Software Foundation. GNU Core Utilities. https://www.gnu.org/software/coreutils/

---

## GNU grep

Purpose:
Pattern matching and text searching.

Main functionalities used:
- Filtering outputs
- Searching log files
- Pattern extraction

Citation:

Free Software Foundation. GNU Grep Manual. https://www.gnu.org/software/grep/

---

## GNU sed

Purpose:
Stream editing and text transformation.

Main functionalities used:
- Text substitution
- Pipeline text manipulation

Citation:

Free Software Foundation. GNU sed Manual. https://www.gnu.org/software/sed/

---

## Bash

Purpose:
Shell scripting and command-line workflow execution.

Main functionalities used:
- Pipeline scripting
- Workflow orchestration
- Command chaining

Citation:

Free Software Foundation. Bash Reference Manual. https://www.gnu.org/software/bash/
