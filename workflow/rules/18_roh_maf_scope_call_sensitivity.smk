# ============================================================
# ROH MAF-scope call sensitivity
# ============================================================
# Final marker-ascertainment sensitivity before production freezing:
# compare the selected ROH parameters on the current joint-MAF panel versus
# the stricter intersection requiring MAF>=0.05 within all six populations.

ROH_MAF_CALL_DIR = "results/roh/maf_scope_call_sensitivity"
ROH_ALLPOP_CALL_PREFIX = ROH_MAF_CALL_DIR + "/all_population_maf_intersection"


rule run_roh_all_population_maf_sensitivity:
    input:
        bed=ROH_COMMON_BED_PREFIX + ".bed",
        bim=ROH_COMMON_BED_PREFIX + ".bim",
        fam=ROH_COMMON_BED_PREFIX + ".fam",
        snplist=ROH_ALLPOP_MAF_SNPLIST
    output:
        hom=ROH_ALLPOP_CALL_PREFIX + ".hom",
        indiv=ROH_ALLPOP_CALL_PREFIX + ".hom.indiv",
        summary=ROH_ALLPOP_CALL_PREFIX + ".hom.summary"
    log:
        "logs/roh/maf_scope_call_sensitivity/all_population_maf_intersection.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_MAF_CALL_DIR} logs/roh/maf_scope_call_sensitivity

        plink \
            --bfile {ROH_COMMON_BED_PREFIX} \
            --extract {input.snplist} \
            --homozyg \
            --homozyg-window-snp 50 \
            --homozyg-snp 50 \
            --homozyg-kb 1500 \
            --homozyg-density 50 \
            --homozyg-gap 500 \
            --homozyg-window-missing 5 \
            --homozyg-window-het 1 \
            --homozyg-window-threshold 0.05 \
            --out {ROH_ALLPOP_CALL_PREFIX} \
            > {log} 2>&1
        """


rule summarize_roh_maf_scope_call_sensitivity:
    input:
        joint_hom=ROH_GAP_HET_DIR + "/gap500_het1.hom",
        allpop_hom=ROH_ALLPOP_CALL_PREFIX + ".hom",
        fam=ROH_COMMON_BED_PREFIX + ".fam",
        annotations=PCA_SAMPLE_ANNOTATIONS,
        script="workflow/scripts/summarize_roh_maf_scope_call_sensitivity.py"
    output:
        individual=ROH_MAF_CALL_DIR + "/individual_summary.tsv",
        population=ROH_MAF_CALL_DIR + "/population_summary.tsv",
        comparison=ROH_MAF_CALL_DIR + "/individual_comparison.tsv",
        correlations=ROH_MAF_CALL_DIR + "/correlations.tsv",
        readme=ROH_MAF_CALL_DIR + "/README.txt"
    params:
        expected_samples=lambda wildcards: config["population_structure"]["roh"]["expected_joint_samples"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --joint-hom {input.joint_hom} \
            --allpop-hom {input.allpop_hom} \
            --fam {input.fam} \
            --annotations {input.annotations} \
            --expected-samples {params.expected_samples} \
            --individual-out {output.individual} \
            --population-out {output.population} \
            --comparison-out {output.comparison} \
            --correlations-out {output.correlations} \
            --readme-out {output.readme}
        """


rule compare_roh_maf_scope_call_sensitivity:
    input:
        ROH_ALLPOP_CALL_PREFIX + ".hom",
        ROH_MAF_CALL_DIR + "/individual_summary.tsv",
        ROH_MAF_CALL_DIR + "/population_summary.tsv",
        ROH_MAF_CALL_DIR + "/individual_comparison.tsv",
        ROH_MAF_CALL_DIR + "/correlations.tsv",
        ROH_MAF_CALL_DIR + "/README.txt"
