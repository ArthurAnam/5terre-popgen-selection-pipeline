# ============================================================
# Cinque Terre Population Genomics and Selection Pipeline
# ============================================================

import os
import yaml

# Version-controlled analytical configuration.
configfile: "config/config.yaml"

# Machine-specific paths are kept outside Git. When present, this file is
# merged on top of config/config.yaml so that rules never hard-code local paths.
LOCAL_CONFIG = "config/config.local.yaml"
if os.path.exists(LOCAL_CONFIG):
    with open(LOCAL_CONFIG, "r") as handle:
        local_config = yaml.safe_load(handle) or {}
    config.update(local_config)


# ============================================================
# Workflow modules
# ============================================================

include: "workflow/rules/00_preqc.smk"
include: "workflow/rules/01_sample_qc.smk"
include: "workflow/rules/02_sample_preqc_chromosome.smk"
include: "workflow/rules/03_sample_preqc_genotype_fields.smk"
include: "workflow/rules/04_sample_preqc_plots.smk"
include: "workflow/rules/05_target_qc_start.smk"
include: "workflow/rules/06_genotype_qc_exploration.smk"
include: "workflow/rules/07_core_qc_filters.smk"
include: "workflow/rules/08_king_relatedness.smk"
include: "workflow/rules/99_provenance.smk"


# ============================================================
# Current reproducible target
# ============================================================
# Individual pre-QC is diagnostic at this stage: no sample is removed
# automatically from the delivered 50-sample dataset.

rule all:
    input:
        "results/preqc/preqc_done.txt",
        "results/sample_preqc/sample_preqc_summary.tsv",
        "results/sample_preqc/sample_chromosome_consistency.tsv",
        "results/sample_preqc/sample_genotype_field_diagnostics.tsv",
        "results/sample_preqc/plots/depth_vs_missingness.png",
        "results/sample_preqc/plots/heterozygosity_vs_allelic_imbalance.png",
        "results/sample_preqc/plots/README.txt",
        "results/qc/00_post_sample_qc_input/post_sample_qc_input_summary.tsv",
        "results/qc/00_post_sample_qc_input/README.txt",
        "results/qc/01_variant_qc_exploration/variant_qc_metric_descriptive_statistics.tsv",
        "results/qc/01_variant_qc_exploration/variant_category_counts_before_filtering.tsv",
        "results/qc/01_variant_qc_exploration/variant_qc_metric_distributions.png",
        "results/qc/01_variant_qc_exploration/README.txt",
        "results/qc/02_genotype_qc_exploration/genotype_qc_metric_descriptive_statistics.tsv",
        "results/qc/02_genotype_qc_exploration/genotype_qc_metric_distributions.png",
        "results/qc/02_genotype_qc_exploration/individual_missingness_before_variant_filtering.tsv",
        "results/qc/02_genotype_qc_exploration/README.txt",
        "results/qc/03_core_filtering/core_variant_filter_step_counts.tsv",
        "results/qc/03_core_filtering/iterative_missingness_qc_history.tsv",
        "results/qc/03_core_filtering/individual_missingness_after_site_filters.tsv",
        "results/qc/03_core_filtering/individuals_excluded_missingness_gt_0.05.txt",
        "results/qc/03_core_filtering/retained_samples_before_hwe.txt",
        "results/qc/03_core_filtering/core_qc_pre_hwe_summary.tsv",
        "results/qc/03_core_filtering/README.txt",
        "results/qc/04_hwe/cinque_terre.qc_filtered.vcf.gz",
        "results/qc/04_hwe/cinque_terre.qc_filtered.vcf.gz.tbi",
        "results/qc/04_hwe/hwe_bonferroni_thresholds.tsv",
        "results/qc/04_hwe/hwe_bonferroni_by_chromosome.tsv",
        "results/qc/04_hwe/hwe_failed_variants.tsv",
        "results/qc/04_hwe/final_qc_summary.tsv",
        "results/qc/04_hwe/README.txt",
        "results/relatedness/king/input/cinque_terre.king.bed",
        "results/relatedness/king/input/cinque_terre.king.bim",
        "results/relatedness/king/input/cinque_terre.king.fam",
        "results/relatedness/king/king_input_dataset_summary.tsv",
        "results/relatedness/king/king.kin0",
        "results/relatedness/king/king_pairwise_kinship.tsv",
        "results/relatedness/king/king_candidate_relatives.tsv",
        "results/relatedness/king/king_relationship_counts.tsv",
        "results/relatedness/king/king_relatedness_summary.tsv",
        "results/relatedness/king/king_individual_kinship_summary.tsv",
        "results/relatedness/king/king_kinship_distribution.png",
        "results/relatedness/king/king_kinship_vs_ibs0.png",
        "results/relatedness/king/king_kinship_vs_ibs0.pdf",
        "results/relatedness/king/king_kinship_heatmap.png",
        "results/relatedness/king/INPUT_README.txt",
        "results/relatedness/king/README.txt",
        "results/provenance/software_versions.tsv"
