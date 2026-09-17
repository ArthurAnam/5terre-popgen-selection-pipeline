# Upstream variant-calling provenance

## Scope

This document records how the Cinque Terre WGS VCF used as input to this repository was generated upstream. Variant calling and VQSR are therefore **provenance**, not steps repeated by the downstream Snakemake workflow.

The distinction between reported upstream information and facts independently verified from the delivered VCF is intentional. It prevents downstream analyses from silently treating historical notes or provider information as newly validated results.

---

## Upstream provenance reported for the dataset

The following description was supplied with the previous analysis material and project notes.

Whole-genome sequencing was performed by Dante Labs using the Illumina NovaSeq X Plus platform. Genomic libraries were prepared with the Illumina DNA PCR-Free Prep Kit, sequencing generated paired-end reads with a reported mean length of approximately 147-148 bp, and mean sequencing coverage was reported to be at least 30x. Raw data were demultiplexed and converted to FASTQ using BCL Convert. Read alignment and initial provider-side processing were performed with the Illumina DRAGEN Bio-IT Platform. Detailed internal DRAGEN parameters were not disclosed.

The downstream joint-calling workflow is reported to have used GATK 4.1.9.0. Individual samples were represented as gVCFs, imported jointly using `GenomicsDBImport`, and genotyped using `GenotypeGVCFs`. The reported reference was GRCh37.p13 with hs37d5 decoy sequences.

Variant Quality Score Recalibration (VQSR) was reported to follow the Broad Institute germline short-variant workflow. Recalibration was performed separately for SNPs and INDELs using standard GATK resource-bundle training/truth resources. For SNPs, these included HapMap, Omni, 1000 Genomes, and dbSNP; INDEL recalibration additionally used Mills and AxiomPoly resources. The reported truth-sensitivity filter level was 99.7%.

The VQSR model used variant-level annotations including QD, MQ, MQRankSum, ReadPosRankSum, FS, SOR, and DP. VQSR-related INFO annotations such as `VQSLOD`, `POSITIVE_TRAIN_SITE`, `NEGATIVE_TRAIN_SITE`, and `culprit` were generated for model training/diagnostic purposes and were not treated as independent biological filters.

After recalibration, the dataset was restricted to biallelic SNPs using GATK `SelectVariants`, excluding multiallelic sites and INDELs. The delivered analysis VCF was then restricted to records passing upstream filtering (`FILTER=PASS`). The person who supplied the VCF explicitly reported that variants failing Variant Recalibration filters had been removed before delivery.

A Ts/Tv ratio of approximately 2.12 was reported in the previous analysis material. This value should be treated as historical provenance until explicitly cross-checked against the current reproducible outputs.

---

## Technical provenance table

| Stage | Software / tool | Version | Reference / resources | Key parameters / thresholds | Input -> output | Provenance status |
|---|---|---:|---|---|---|---|
| Sequencing & base calling | Illumina NovaSeq X Plus | not disclosed | n/a | Reported >=30x mean coverage; PE reads ~147-148 bp | Blood/DBS -> sequence data | Provider/project notes |
| Demultiplexing / FASTQ generation | BCL Convert | not disclosed | n/a | Provider standard protocol | BCL -> FASTQ | Provider report |
| Alignment / provider-side processing | DRAGEN Bio-IT Platform | not disclosed | Reported hs37d5 | Internal parameters not disclosed | FASTQ -> aligned data / provider intermediates | Provider report |
| gVCF joint import | GATK `GenomicsDBImport` | 4.1.9.0 | Reported hs37d5 | Sample map; per-chromosome workspace | Per-sample gVCFs -> GenomicsDB | Previous workflow / VCF-header provenance |
| Joint genotyping | GATK `GenotypeGVCFs` | 4.1.9.0 | Reported hs37d5 | Workflow defaults / historical configuration | GenomicsDB -> multisample VCF | Previous workflow provenance |
| SNP extraction | GATK `SelectVariants` | 4.1.9.0 | same reference | `--select-type-to-include SNP`; `--restrict-alleles-to BIALLELIC` | Joint VCF -> SNP-only VCF | Previous workflow / header provenance |
| SNP VQSR application | GATK `ApplyVQSR` | 4.1.9.0 | HapMap, Omni, 1000G, dbSNP | mode=SNP; truth-sensitivity-filter-level 99.7 | Raw SNP VCF -> recalibrated SNP VCF | Previous workflow provenance |
| Delivered input definition | PASS-only biallelic SNPs | n/a | n/a | retain `FILTER=PASS` | upstream VCF -> `all_samples_snp.vcf.gz` | Historical report + current pre-QC verification |

---

## Facts independently verified by the current pipeline

The reproducible pre-QC rule in `workflow/rules/00_preqc.smk` has now validated the delivered input `all_samples_snp.vcf.gz` without modifying it.

Current verified properties:

- 50 input samples
- 11,490,646 VCF records
- 11,490,646 SNPs
- 0 INDELs
- 0 MNPs
- 0 other variant classes
- 0 multiallelic sites
- 0 multiallelic SNP sites
- 11,490,646 records with `FILTER=PASS`
- 0 non-PASS records
- indexed contigs: chromosomes `1`-`22`, `X`, and `Y`

These observations establish the reproducible starting point for the downstream QC workflow. In particular, the pipeline does **not** repeat VQSR or SNP/biallelic filtering merely because those operations occurred upstream; instead, it verifies the delivered file and records its properties.

The downstream population-genomics workflow will subsequently restrict analyses to autosomes where required.

---

## Sample-count clarification

The delivered input contains **50 samples**. Earlier project documentation referred to **46 samples**, but that appears to represent a historical post-QC dataset rather than the original delivered input.

The current pipeline therefore starts from 50 individuals and does not pre-impose a final sample count. The historical reduction from 50 to 46 will be reconstructed from reproducible sample-QC evidence before any four individuals are excluded.

---

## Reproducibility boundary

This repository begins from the delivered post-VQSR SNP VCF. Upstream sequencing, alignment, gVCF generation, joint genotyping, VQSR model fitting/application, and initial biallelic-SNP definition are documented here for provenance but are not rerun by default.

The reproducible downstream boundary is therefore:

`all_samples_snp.vcf.gz` (50 samples, 11,490,646 PASS SNPs) -> pre-QC verification -> target sample/variant QC -> population-genomic analyses -> phasing and selection analyses.

Where an upstream parameter is unavailable or supported only by provider/project notes, it remains explicitly labelled as such rather than inferred.
