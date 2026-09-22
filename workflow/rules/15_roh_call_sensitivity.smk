# ROH call sensitivity: population-history vs light-VIF

ROH_CALL_SENS_DIR = "results/roh/call_sensitivity"
ROH_PRIMARY_PREFIX = ROH_CALL_SENS_DIR + "/population_history_unpruned"
ROH_VIF_PREFIX = ROH_CALL_SENS_DIR + "/light_vif_howrigan"

rule run_roh_population_history_candidate:
    input:
        bed=ROH_COMMON_BED_PREFIX + ".bed",
        bim=ROH_COMMON_BED_PREFIX + ".bim",
        fam=ROH_COMMON_BED_PREFIX + ".fam"
    output:
        hom=ROH_PRIMARY_PREFIX + ".hom",
        indiv=ROH_PRIMARY_PREFIX + ".hom.indiv",
        summary=ROH_PRIMARY_PREFIX + ".hom.summary"
    log:
        "logs/roh/call_sensitivity/population_history_unpruned.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_CALL_SENS_DIR} logs/roh/call_sensitivity
        plink \
            --bfile {ROH_COMMON_BED_PREFIX} \
            --homozyg \
            --homozyg-window-snp 50 \
            --homozyg-snp 50 \
            --homozyg-kb 1500 \
            --homozyg-density 50 \
            --homozyg-gap 1000 \
            --homozyg-window-missing 5 \
            --homozyg-window-het 1 \
            --homozyg-window-threshold 0.05 \
            --out {ROH_PRIMARY_PREFIX} \
            > {log} 2>&1
        """

rule run_roh_light_vif_sensitivity:
    input:
        bed=ROH_COMMON_BED_PREFIX + ".bed",
        bim=ROH_COMMON_BED_PREFIX + ".bim",
        fam=ROH_COMMON_BED_PREFIX + ".fam",
        prune_in=ROH_LIGHT_VIF_PREFIX + ".prune.in"
    output:
        hom=ROH_VIF_PREFIX + ".hom",
        indiv=ROH_VIF_PREFIX + ".hom.indiv",
        summary=ROH_VIF_PREFIX + ".hom.summary"
    log:
        "logs/roh/call_sensitivity/light_vif_howrigan.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_CALL_SENS_DIR} logs/roh/call_sensitivity
        plink \
            --bfile {ROH_COMMON_BED_PREFIX} \
            --extract {input.prune_in} \
            --homozyg \
            --homozyg-window-snp 65 \
            --homozyg-snp 65 \
            --homozyg-kb 10 \
            --homozyg-density 200 \
            --homozyg-gap 500 \
            --homozyg-window-missing 3 \
            --homozyg-window-het 0 \
            --homozyg-window-threshold 0.05 \
            --out {ROH_VIF_PREFIX} \
            > {log} 2>&1
        """

rule summarize_roh_call_sensitivity:
    input:
        primary_hom=ROH_PRIMARY_PREFIX + ".hom",
        vif_hom=ROH_VIF_PREFIX + ".hom",
        fam=ROH_COMMON_BED_PREFIX + ".fam",
        annotations=PCA_SAMPLE_ANNOTATIONS,
        script="workflow/scripts/summarize_roh_call_sensitivity.py"
    output:
        segments=ROH_CALL_SENS_DIR + "/long_roh_segments.tsv",
        individual=ROH_CALL_SENS_DIR + "/individual_long_roh_summary.tsv",
        population=ROH_CALL_SENS_DIR + "/population_long_roh_summary.tsv",
        comparison=ROH_CALL_SENS_DIR + "/individual_framework_comparison.tsv",
        readme=ROH_CALL_SENS_DIR + "/README.txt"
    params:
        min_kb=1500,
        expected_samples=lambda wildcards: config["population_structure"]["roh"]["expected_joint_samples"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --primary-hom {input.primary_hom} \
            --vif-hom {input.vif_hom} \
            --fam {input.fam} \
            --annotations {input.annotations} \
            --min-kb {params.min_kb} \
            --expected-samples {params.expected_samples} \
            --segments-out {output.segments} \
            --individual-out {output.individual} \
            --population-out {output.population} \
            --comparison-out {output.comparison} \
            --readme-out {output.readme}
        """

rule compare_roh_call_sensitivity:
    input:
        ROH_PRIMARY_PREFIX + ".hom",
        ROH_PRIMARY_PREFIX + ".hom.indiv",
        ROH_VIF_PREFIX + ".hom",
        ROH_VIF_PREFIX + ".hom.indiv",
        ROH_CALL_SENS_DIR + "/long_roh_segments.tsv",
        ROH_CALL_SENS_DIR + "/individual_long_roh_summary.tsv",
        ROH_CALL_SENS_DIR + "/population_long_roh_summary.tsv",
        ROH_CALL_SENS_DIR + "/individual_framework_comparison.tsv",
        ROH_CALL_SENS_DIR + "/README.txt"
