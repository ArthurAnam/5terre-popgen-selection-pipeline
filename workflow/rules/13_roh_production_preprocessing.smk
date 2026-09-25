# ============================================================
# Production ROH preprocessing
# ============================================================

ROH_PREPROCESSING_DIR = "results/roh/preprocessing"
ROH_COMMON_BED_PREFIX = ROH_PREPROCESSING_DIR + "/ct_1kg_eur.maf0.05"


rule prepare_roh_production_panel:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam"
    output:
        bed=ROH_COMMON_BED_PREFIX + ".bed",
        bim=ROH_COMMON_BED_PREFIX + ".bim",
        fam=ROH_COMMON_BED_PREFIX + ".fam",
        summary=ROH_PREPROCESSING_DIR + "/production_panel_summary.tsv"
    params:
        maf=lambda wc: config["population_structure"]["roh"]["maf_threshold"],
        expected_variants=lambda wc: config["population_structure"]["roh"]["expected_joint_maf0.05_variants"],
        expected_samples=lambda wc: config["population_structure"]["roh"]["expected_joint_samples"]
    log:
        "logs/roh/preprocessing/production_panel.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_PREPROCESSING_DIR} logs/roh/preprocessing

        tmp_prefix="{ROH_PREPROCESSING_DIR}/.production_panel.$$"
        trap 'rm -f "$tmp_prefix".*' EXIT

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --maf {params.maf} \
            --make-bed \
            --out "$tmp_prefix" \
            > {log} 2>&1

        n_variants=$(wc -l < "$tmp_prefix.bim")
        n_samples=$(wc -l < "$tmp_prefix.fam")

        [ "$n_variants" -eq "{params.expected_variants}" ] || {
            echo "ERROR: expected {params.expected_variants} ROH markers, found $n_variants" >&2
            exit 1
        }
        [ "$n_samples" -eq "{params.expected_samples}" ] || {
            echo "ERROR: expected {params.expected_samples} ROH samples, found $n_samples" >&2
            exit 1
        }

        mv "$tmp_prefix.bed" {output.bed}
        mv "$tmp_prefix.bim" {output.bim}
        mv "$tmp_prefix.fam" {output.fam}

        {
            printf 'metric\tvalue\n'
            printf 'maf_threshold\t%s\n' "{params.maf}"
            printf 'variants\t%s\n' "$n_variants"
            printf 'samples\t%s\n' "$n_samples"
            printf 'ld_pruning\tNO\n'
        } > {output.summary}
        """
