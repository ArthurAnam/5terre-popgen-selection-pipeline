# ============================================================
# ROH/FROH population distributions and CT-vs-reference contrasts
# ============================================================

ROH_DIST_DIR = "results/roh/population_distributions"

rule analyze_roh_population_distributions:
    input:
        individual=ROH_PROD_DIR + "/individual_roh_froh.tsv",
        script="workflow/scripts/analyze_roh_population_distributions.py"
    output:
        descriptive=ROH_DIST_DIR + "/population_distribution_summary.tsv",
        tests=ROH_DIST_DIR + "/ct_vs_reference_tests.tsv",
        froh_png=ROH_DIST_DIR + "/froh_ge1.5mb_by_population.png",
        froh_pdf=ROH_DIST_DIR + "/froh_ge1.5mb_by_population.pdf",
        long_png=ROH_DIST_DIR + "/froh_ge5mb_by_population.png",
        long_pdf=ROH_DIST_DIR + "/froh_ge5mb_by_population.pdf",
        scatter_png=ROH_DIST_DIR + "/nroh_vs_total_roh.png",
        scatter_pdf=ROH_DIST_DIR + "/nroh_vs_total_roh.pdf",
        readme=ROH_DIST_DIR + "/README.txt"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_DIST_DIR}
        python {input.script} \
            --individual {input.individual} \
            --descriptive-out {output.descriptive} \
            --tests-out {output.tests} \
            --froh-png {output.froh_png} \
            --froh-pdf {output.froh_pdf} \
            --long-png {output.long_png} \
            --long-pdf {output.long_pdf} \
            --scatter-png {output.scatter_png} \
            --scatter-pdf {output.scatter_pdf} \
            --readme-out {output.readme}
        """

rule summarize_roh_population_distributions:
    input:
        ROH_DIST_DIR + "/population_distribution_summary.tsv",
        ROH_DIST_DIR + "/ct_vs_reference_tests.tsv",
        ROH_DIST_DIR + "/froh_ge1.5mb_by_population.png",
        ROH_DIST_DIR + "/froh_ge1.5mb_by_population.pdf",
        ROH_DIST_DIR + "/froh_ge5mb_by_population.png",
        ROH_DIST_DIR + "/froh_ge5mb_by_population.pdf",
        ROH_DIST_DIR + "/nroh_vs_total_roh.png",
        ROH_DIST_DIR + "/nroh_vs_total_roh.pdf",
        ROH_DIST_DIR + "/README.txt"
