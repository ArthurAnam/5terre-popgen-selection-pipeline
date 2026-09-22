# ============================================================
# ROH preprocessing audit: common-marker density before calling ROH
# ============================================================

ROH_PREPROCESSING_DIR = "results/roh/preprocessing"
ROH_MAF_SNPLIST = ROH_PREPROCESSING_DIR + "/ct_1kg_eur.roh.maf0.05.snplist"


rule roh_joint_maf_variant_list:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam"
    output:
        ROH_MAF_SNPLIST
    params:
        maf=lambda wildcards: config["population_structure"]["roh"]["maf_threshold"]
    log:
        "logs/roh/preprocessing/maf_filter.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_PREPROCESSING_DIR} logs/roh/preprocessing

        tmp_prefix="{ROH_PREPROCESSING_DIR}/.maf_filter.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --maf {params.maf} \
            --write-snplist \
            --out "$tmp_prefix" \
            > {log} 2>&1

        mv "$tmp_prefix.snplist" {output}
        """


rule audit_roh_marker_density:
    input:
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        snplist=ROH_MAF_SNPLIST,
        script="workflow/scripts/summarize_roh_marker_density.py"
    output:
        by_chr=ROH_PREPROCESSING_DIR + "/marker_density_by_chromosome.tsv",
        gaps=ROH_PREPROCESSING_DIR + "/intermarker_gap_thresholds.tsv",
        summary=ROH_PREPROCESSING_DIR + "/roh_marker_density_summary.tsv",
        readme=ROH_PREPROCESSING_DIR + "/README.txt"
    params:
        maf=lambda wildcards: config["population_structure"]["roh"]["maf_threshold"],
        expected_variants=lambda wildcards: config["population_structure"]["roh"]["expected_joint_maf0.05_variants"],
        expected_samples=lambda wildcards: config["population_structure"]["roh"]["expected_joint_samples"],
        gap_thresholds=lambda wildcards: ",".join(
            str(x) for x in config["population_structure"]["roh"]["marker_density_audit"]["gap_thresholds_kb"]
        )
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --pvar {input.pvar} \
            --snplist {input.snplist} \
            --maf-threshold {params.maf} \
            --expected-variants {params.expected_variants} \
            --expected-samples {params.expected_samples} \
            --gap-thresholds-kb {params.gap_thresholds} \
            --by-chromosome-out {output.by_chr} \
            --gap-thresholds-out {output.gaps} \
            --summary-out {output.summary} \
            --readme-out {output.readme}
        """


rule audit_roh_preprocessing:
    input:
        ROH_MAF_SNPLIST,
        ROH_PREPROCESSING_DIR + "/marker_density_by_chromosome.tsv",
        ROH_PREPROCESSING_DIR + "/intermarker_gap_thresholds.tsv",
        ROH_PREPROCESSING_DIR + "/roh_marker_density_summary.tsv",
        ROH_PREPROCESSING_DIR + "/README.txt"
