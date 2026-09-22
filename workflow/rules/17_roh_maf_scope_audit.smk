# ============================================================
# ROH MAF-scope audit
# ============================================================
# Compare the current joint MAF>=0.05 common-marker panel with a deliberately
# stricter sensitivity panel requiring MAF>=0.05 within every population.
# This stage audits marker count/density only; it does not yet replace the
# primary ROH panel or launch an alternate production ROH call.

ROH_MAF_SCOPE_DIR = "results/roh/maf_scope_sensitivity"
ROH_MAF_SCOPE_POPS = ["CT", "CEU", "FIN", "GBR", "IBS", "TSI"]
ROH_ALLPOP_MAF_SNPLIST = ROH_MAF_SCOPE_DIR + "/all_populations.maf0.05.intersection.snplist"


rule roh_maf_scope_population_keep_files:
    input:
        fam=ROH_COMMON_BED_PREFIX + ".fam",
        annotations=PCA_SAMPLE_ANNOTATIONS,
        script="workflow/scripts/make_roh_population_keep_files.py"
    output:
        expand(ROH_MAF_SCOPE_DIR + "/keep/{population}.keep", population=ROH_MAF_SCOPE_POPS)
    params:
        outdir=ROH_MAF_SCOPE_DIR + "/keep",
        expected="CT:46,CEU:99,FIN:99,GBR:91,IBS:107,TSI:107"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --fam {input.fam} \
            --annotations {input.annotations} \
            --out-dir {params.outdir} \
            --expected-counts {params.expected}
        """


rule roh_population_maf_snplist:
    input:
        bed=ROH_COMMON_BED_PREFIX + ".bed",
        bim=ROH_COMMON_BED_PREFIX + ".bim",
        fam=ROH_COMMON_BED_PREFIX + ".fam",
        keep=ROH_MAF_SCOPE_DIR + "/keep/{population}.keep"
    output:
        ROH_MAF_SCOPE_DIR + "/maf_lists/{population}.maf0.05.snplist"
    log:
        "logs/roh/maf_scope/{population}.maf0.05.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {ROH_MAF_SCOPE_DIR}/maf_lists logs/roh/maf_scope
        tmp="{ROH_MAF_SCOPE_DIR}/maf_lists/.{wildcards.population}.$$"
        trap 'rm -f "$tmp".*' EXIT

        plink2 \
            --threads {threads} \
            --bfile {ROH_COMMON_BED_PREFIX} \
            --keep {input.keep} \
            --maf 0.05 \
            --write-snplist \
            --out "$tmp" \
            > {log} 2>&1

        mv "$tmp.snplist" {output}
        """


rule roh_all_population_maf_intersection:
    input:
        joint=ROH_MAF_SNPLIST,
        lists=expand(
            ROH_MAF_SCOPE_DIR + "/maf_lists/{population}.maf0.05.snplist",
            population=ROH_MAF_SCOPE_POPS
        ),
        script="workflow/scripts/intersect_roh_population_maf_snplists.py"
    output:
        snplist=ROH_ALLPOP_MAF_SNPLIST,
        counts=ROH_MAF_SCOPE_DIR + "/per_population_maf_counts.tsv",
        summary=ROH_MAF_SCOPE_DIR + "/intersection_summary.tsv"
    params:
        pop_args=" ".join(
            f"--population-list {population}={ROH_MAF_SCOPE_DIR}/maf_lists/{population}.maf0.05.snplist"
            for population in ROH_MAF_SCOPE_POPS
        )
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --joint-snplist {input.joint} \
            {params.pop_args} \
            --intersection-out {output.snplist} \
            --counts-out {output.counts} \
            --summary-out {output.summary}
        """


rule audit_roh_all_population_maf_density:
    input:
        pvar=PCA_HARMONIZED_PREFIX + ".pvar",
        snplist=ROH_ALLPOP_MAF_SNPLIST,
        script="workflow/scripts/summarize_roh_marker_density.py"
    output:
        by_chr=ROH_MAF_SCOPE_DIR + "/marker_density_by_chromosome.tsv",
        gaps=ROH_MAF_SCOPE_DIR + "/intermarker_gap_thresholds.tsv",
        summary=ROH_MAF_SCOPE_DIR + "/marker_density_summary.tsv",
        readme=ROH_MAF_SCOPE_DIR + "/README.txt"
    params:
        expected_samples=lambda wildcards: config["population_structure"]["roh"]["expected_joint_samples"],
        gap_thresholds=lambda wildcards: ",".join(
            str(x) for x in config["population_structure"]["roh"]["marker_density_audit"]["gap_thresholds_kb"]
        )
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        n_variants=$(wc -l < {input.snplist})

        python {input.script} \
            --pvar {input.pvar} \
            --snplist {input.snplist} \
            --maf-threshold 0.05 \
            --expected-variants "$n_variants" \
            --expected-samples {params.expected_samples} \
            --panel-label all_populations_individually_maf_ge_0.05_intersection \
            --ld-pruning-label NO \
            --gap-thresholds-kb {params.gap_thresholds} \
            --by-chromosome-out {output.by_chr} \
            --gap-thresholds-out {output.gaps} \
            --summary-out {output.summary} \
            --readme-out {output.readme}
        """


rule audit_roh_maf_scope_sensitivity:
    input:
        ROH_ALLPOP_MAF_SNPLIST,
        ROH_MAF_SCOPE_DIR + "/per_population_maf_counts.tsv",
        ROH_MAF_SCOPE_DIR + "/intersection_summary.tsv",
        ROH_MAF_SCOPE_DIR + "/marker_density_by_chromosome.tsv",
        ROH_MAF_SCOPE_DIR + "/intermarker_gap_thresholds.tsv",
        ROH_MAF_SCOPE_DIR + "/marker_density_summary.tsv",
        ROH_MAF_SCOPE_DIR + "/README.txt"
