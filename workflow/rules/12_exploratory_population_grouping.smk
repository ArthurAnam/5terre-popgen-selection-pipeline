# ============================================================
# Exploratory population grouping: pairwise Hudson FST
# ============================================================
#
# Explicit target only. This branch does not change any population label used
# elsewhere in the workflow. It reuses the joint common-variant panel after
# the predefined high-LD mask, before LD pruning.

GROUPING_DIR = "results/population_structure/exploratory_grouping"
GROUPING_WITHIN = GROUPING_DIR + "/population_clusters.txt"
GROUPING_FST_PREFIX = GROUPING_DIR + "/pairwise_hudson"
GROUPING_FST_SUMMARY = GROUPING_FST_PREFIX + ".fst.summary"


rule make_exploratory_grouping_clusters:
    input:
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        annotations=PCA_SAMPLE_ANNOTATIONS,
        script="workflow/scripts/make_plink_within_from_annotations.py"
    output:
        GROUPING_WITHIN
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {GROUPING_DIR}
        python {input.script} \
            --psam {input.psam} \
            --annotations {input.annotations} \
            --output {output}
        """


rule exploratory_pairwise_hudson_fst:
    input:
        pgen=PCA_HARMONIZED_PREFIX + ".pgen",
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        psam=PCA_HARMONIZED_PREFIX + ".psam",
        variants=PCA_HIGHLD_SNPLIST,
        clusters=GROUPING_WITHIN
    output:
        GROUPING_FST_SUMMARY
    log:
        "logs/population_structure/exploratory_grouping/pairwise_hudson_fst.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {GROUPING_DIR} logs/population_structure/exploratory_grouping

        plink2 \
            --threads {threads} \
            --pfile {PCA_HARMONIZED_PREFIX} \
            --extract {input.variants} \
            --within {input.clusters} PCA_POP \
            --fst PCA_POP method=hudson \
            --out {GROUPING_FST_PREFIX} \
            > {log} 2>&1
        """


rule summarize_exploratory_population_grouping:
    input:
        fst=GROUPING_FST_SUMMARY,
        script="workflow/scripts/summarize_population_grouping.py"
    output:
        matrix=GROUPING_DIR + "/pairwise_hudson_fst_matrix.tsv",
        ranked=GROUPING_DIR + "/pairwise_hudson_fst_ranked.tsv",
        merges=GROUPING_DIR + "/average_linkage_merges.tsv",
        dendrogram=GROUPING_DIR + "/average_linkage_dendrogram.png",
        readme=GROUPING_DIR + "/README.txt"
    params:
        populations=lambda wildcards: ",".join(
            config["population_structure"]["exploratory_population_grouping"]["populations"]
        )
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --fst-summary {input.fst} \
            --populations {params.populations} \
            --matrix-out {output.matrix} \
            --ranked-out {output.ranked} \
            --merges-out {output.merges} \
            --dendrogram-out {output.dendrogram} \
            --readme-out {output.readme}
        """


rule explore_population_grouping:
    input:
        GROUPING_FST_SUMMARY,
        GROUPING_DIR + "/pairwise_hudson_fst_matrix.tsv",
        GROUPING_DIR + "/pairwise_hudson_fst_ranked.tsv",
        GROUPING_DIR + "/average_linkage_merges.tsv",
        GROUPING_DIR + "/average_linkage_dendrogram.png",
        GROUPING_DIR + "/README.txt"
