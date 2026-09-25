# ============================================================
# Production ROH and FROH
# ============================================================

ROH_PROD_DIR = "results/roh/production"
ROH_PROD_PREFIX = ROH_PROD_DIR + "/ct_1kg_eur.primary_roh"

rule run_roh_production:
    input:
        bed=ROH_COMMON_BED_PREFIX + ".bed",
        bim=ROH_COMMON_BED_PREFIX + ".bim",
        fam=ROH_COMMON_BED_PREFIX + ".fam"
    output:
        hom=ROH_PROD_PREFIX + ".hom",
        indiv=ROH_PROD_PREFIX + ".hom.indiv",
        summary=ROH_PROD_PREFIX + ".hom.summary"
    params:
        window_snp=lambda wc: config["population_structure"]["roh"]["production"]["homozyg_window_snp"],
        min_snp=lambda wc: config["population_structure"]["roh"]["production"]["homozyg_snp"],
        min_kb=lambda wc: config["population_structure"]["roh"]["production"]["homozyg_kb"],
        density=lambda wc: config["population_structure"]["roh"]["production"]["homozyg_density"],
        gap=lambda wc: config["population_structure"]["roh"]["production"]["homozyg_gap"],
        window_missing=lambda wc: config["population_structure"]["roh"]["production"]["homozyg_window_missing"],
        window_het=lambda wc: config["population_structure"]["roh"]["production"]["homozyg_window_het"],
        window_threshold=lambda wc: config["population_structure"]["roh"]["production"]["homozyg_window_threshold"]
    log:
        "logs/roh/production/primary_roh.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_PROD_DIR} logs/roh/production

        plink \
            --bfile {ROH_COMMON_BED_PREFIX} \
            --homozyg \
            --homozyg-window-snp {params.window_snp} \
            --homozyg-snp {params.min_snp} \
            --homozyg-kb {params.min_kb} \
            --homozyg-density {params.density} \
            --homozyg-gap {params.gap} \
            --homozyg-window-missing {params.window_missing} \
            --homozyg-window-het {params.window_het} \
            --homozyg-window-threshold {params.window_threshold} \
            --out {ROH_PROD_PREFIX} \
            > {log} 2>&1
        """

rule summarize_roh_production:
    input:
        hom=ROH_PROD_PREFIX + ".hom",
        fam=ROH_COMMON_BED_PREFIX + ".fam",
        annotations=PCA_SAMPLE_ANNOTATIONS,
        script="workflow/scripts/summarize_roh_production.py"
    output:
        segments=ROH_PROD_DIR + "/roh_segments.tsv",
        individual=ROH_PROD_DIR + "/individual_roh_froh.tsv",
        population=ROH_PROD_DIR + "/population_roh_froh_summary.tsv",
        readme=ROH_PROD_DIR + "/README.txt"
    params:
        denominator_bp=lambda wc: config["population_structure"]["roh"]["production"]["froh_denominator_bp"],
        expected_samples=lambda wc: config["population_structure"]["roh"]["expected_joint_samples"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --hom {input.hom} \
            --fam {input.fam} \
            --annotations {input.annotations} \
            --denominator-bp {params.denominator_bp} \
            --expected-samples {params.expected_samples} \
            --segments-out {output.segments} \
            --individual-out {output.individual} \
            --population-out {output.population} \
            --readme-out {output.readme}
        """

rule produce_roh_froh:
    input:
        ROH_PROD_PREFIX + ".hom",
        ROH_PROD_PREFIX + ".hom.indiv",
        ROH_PROD_PREFIX + ".hom.summary",
        ROH_PROD_DIR + "/roh_segments.tsv",
        ROH_PROD_DIR + "/individual_roh_froh.tsv",
        ROH_PROD_DIR + "/population_roh_froh_summary.tsv",
        ROH_PROD_DIR + "/README.txt"
