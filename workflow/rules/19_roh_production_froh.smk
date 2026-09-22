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
            --homozyg-window-snp 50 \
            --homozyg-snp 50 \
            --homozyg-kb 1500 \
            --homozyg-density 50 \
            --homozyg-gap 500 \
            --homozyg-window-missing 5 \
            --homozyg-window-het 1 \
            --homozyg-window-threshold 0.05 \
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
        denominator_bp=2770000000,
        expected_samples=lambda wildcards: config["population_structure"]["roh"]["expected_joint_samples"]
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
