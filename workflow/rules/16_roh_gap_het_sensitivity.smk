# ============================================================
# ROH gap x heterozygote sensitivity on the unpruned primary panel
# ============================================================
# Hold marker panel, MAF, physical threshold, SNP count, density and window
# threshold fixed. Vary only maximum internal gap and tolerated heterozygotes
# so their effects can be attributed directly.

ROH_GAP_HET_DIR = "results/roh/gap_het_sensitivity"
ROH_GAP_VALUES = ["100", "500", "1000"]
ROH_HET_VALUES = ["0", "1"]

wildcard_constraints:
    roh_gap="100|500|1000",
    roh_het="0|1"


rule run_roh_gap_het_sensitivity:
    input:
        bed=ROH_COMMON_BED_PREFIX + ".bed",
        bim=ROH_COMMON_BED_PREFIX + ".bim",
        fam=ROH_COMMON_BED_PREFIX + ".fam"
    output:
        hom=ROH_GAP_HET_DIR + "/gap{roh_gap}_het{roh_het}.hom",
        indiv=ROH_GAP_HET_DIR + "/gap{roh_gap}_het{roh_het}.hom.indiv",
        summary=ROH_GAP_HET_DIR + "/gap{roh_gap}_het{roh_het}.hom.summary"
    log:
        "logs/roh/gap_het_sensitivity/gap{roh_gap}_het{roh_het}.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_GAP_HET_DIR} logs/roh/gap_het_sensitivity

        plink \
            --bfile {ROH_COMMON_BED_PREFIX} \
            --homozyg \
            --homozyg-window-snp 50 \
            --homozyg-snp 50 \
            --homozyg-kb 1500 \
            --homozyg-density 50 \
            --homozyg-gap {wildcards.roh_gap} \
            --homozyg-window-missing 5 \
            --homozyg-window-het {wildcards.roh_het} \
            --homozyg-window-threshold 0.05 \
            --out {ROH_GAP_HET_DIR}/gap{wildcards.roh_gap}_het{wildcards.roh_het} \
            > {log} 2>&1
        """


rule summarize_roh_gap_het_sensitivity:
    input:
        hom=expand(
            ROH_GAP_HET_DIR + "/gap{roh_gap}_het{roh_het}.hom",
            roh_gap=ROH_GAP_VALUES,
            roh_het=ROH_HET_VALUES
        ),
        fam=ROH_COMMON_BED_PREFIX + ".fam",
        annotations=PCA_SAMPLE_ANNOTATIONS,
        script="workflow/scripts/summarize_roh_gap_het_sensitivity.py"
    output:
        individual=ROH_GAP_HET_DIR + "/individual_summary.tsv",
        population=ROH_GAP_HET_DIR + "/population_summary.tsv",
        comparison=ROH_GAP_HET_DIR + "/comparison_to_gap1000_het1.tsv",
        correlations=ROH_GAP_HET_DIR + "/correlations_to_gap1000_het1.tsv",
        readme=ROH_GAP_HET_DIR + "/README.txt"
    params:
        expected_samples=lambda wildcards: config["population_structure"]["roh"]["expected_joint_samples"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --input-dir {ROH_GAP_HET_DIR} \
            --fam {input.fam} \
            --annotations {input.annotations} \
            --expected-samples {params.expected_samples} \
            --individual-out {output.individual} \
            --population-out {output.population} \
            --comparison-out {output.comparison} \
            --correlations-out {output.correlations} \
            --readme-out {output.readme}
        """


rule audit_roh_gap_het_sensitivity:
    input:
        expand(
            ROH_GAP_HET_DIR + "/gap{roh_gap}_het{roh_het}.hom",
            roh_gap=ROH_GAP_VALUES,
            roh_het=ROH_HET_VALUES
        ),
        ROH_GAP_HET_DIR + "/individual_summary.tsv",
        ROH_GAP_HET_DIR + "/population_summary.tsv",
        ROH_GAP_HET_DIR + "/comparison_to_gap1000_het1.tsv",
        ROH_GAP_HET_DIR + "/correlations_to_gap1000_het1.tsv",
        ROH_GAP_HET_DIR + "/README.txt"
