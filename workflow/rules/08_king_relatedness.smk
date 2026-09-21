# ============================================================
# KING relatedness analysis: Cinque Terre only
# ============================================================
#
# Architecture:
# 1. create a persistent PLINK BED/BIM/FAM dataset once from the stable pre-HWE
#    target QC dataset;
# 2. run the native KING estimator and keep its raw .kin0 output;
# 3. summarize/plot in a separate lightweight rule;
# 4. generate explanatory README files in documentation-only rules.
#
# This separation avoids rebuilding millions of genotypes when only plots,
# interpretation code, or documentation change.


rule prepare_king_relatedness_dataset:
    input:
        vcf="results/qc/03_core_filtering/cinque_terre.core_qc.pre_hwe.vcf.gz",
        index="results/qc/03_core_filtering/cinque_terre.core_qc.pre_hwe.vcf.gz.tbi"
    output:
        bed="results/relatedness/king/input/cinque_terre.king.bed",
        bim="results/relatedness/king/input/cinque_terre.king.bim",
        fam="results/relatedness/king/input/cinque_terre.king.fam"
    log:
        "logs/relatedness/king/prepare_king_input.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/relatedness/king/input logs/relatedness/king

        plink2 \
            --threads {threads} \
            --vcf {input.vcf} \
            --double-id \
            --set-all-var-ids '@:#:$r:$a' \
            --make-bed \
            --out results/relatedness/king/input/cinque_terre.king \
            > {log} 2>&1
        """


rule summarize_king_input_dataset:
    input:
        bim="results/relatedness/king/input/cinque_terre.king.bim",
        fam="results/relatedness/king/input/cinque_terre.king.fam"
    output:
        "results/relatedness/king/king_input_dataset_summary.tsv"
    shell:
        r"""
        set -euo pipefail
        n_samples=$(wc -l < {input.fam})
        n_variants=$(wc -l < {input.bim})

        {{
            printf 'metric\tvalue\n'
            printf 'source_dataset\tstable_pre_HWE_Cinque_Terre_QC\n'
            printf 'samples\t%s\n' "$n_samples"
            printf 'variants\t%s\n' "$n_variants"
            printf 'autosomes_only\tYES\n'
            printf 'biallelic_SNPs_only\tYES\n'
            printf 'MAF_filter_for_KING\tNO\n'
            printf 'LD_pruning_for_KING\tNO\n'
            printf 'family_ID_assignment\tunique_FID_per_sample_via_PLINK2_double_id\n'
        }} > {output}
        """


rule write_king_input_readme:
    input:
        summary="results/relatedness/king/king_input_dataset_summary.tsv"
    output:
        "results/relatedness/king/INPUT_README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
KING input dataset
==================

The KING dataset is derived from the stable Cinque Terre QC dataset before HWE
filtering.

No extra population-genetic filter is added specifically for KING:
- no MAF filter;
- no LD pruning.

HWE filtering is intentionally not used to define the KING input. Relatedness
is evaluated as a sample-level QC/property of the retained cohort, while KING
recommends retaining genome-wide SNPs that already passed QC rather than
pruning or filtering otherwise good SNPs.

PLINK2 converts the stable pre-HWE QC VCF to PLINK BED/BIM/FAM format. Each
individual is assigned a unique family ID (--double-id), so all 46 choose 2
pairwise comparisons are treated as between-family comparisons and are expected
in KING's .kin0 output.

BED/BIM/FAM are persistent branch-specific inputs. Keeping them prevents
reconstruction of the large pre-HWE VCF when only KING analyses, plots, or
documentation are changed.
EOF
        """


rule run_king_kinship:
    input:
        bed="results/relatedness/king/input/cinque_terre.king.bed",
        bim="results/relatedness/king/input/cinque_terre.king.bim",
        fam="results/relatedness/king/input/cinque_terre.king.fam"
    output:
        kin0="results/relatedness/king/king.kin0"
    log:
        "logs/relatedness/king/king.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/relatedness/king logs/relatedness/king

        king \
            -b {input.bed} \
            --kinship \
            --cpus {threads} \
            --prefix results/relatedness/king/king \
            > {log} 2>&1

        test -s {output.kin0}
        """


rule summarize_and_plot_king_relatedness:
    input:
        kin0="results/relatedness/king/king.kin0",
        fam="results/relatedness/king/input/cinque_terre.king.fam",
        script="workflow/scripts/plot_king_relatedness.py"
    output:
        pairs="results/relatedness/king/king_pairwise_kinship.tsv",
        candidates="results/relatedness/king/king_candidate_relatives.tsv",
        counts="results/relatedness/king/king_relationship_counts.tsv",
        summary="results/relatedness/king/king_relatedness_summary.tsv",
        individual_summary="results/relatedness/king/king_individual_kinship_summary.tsv",
        hist="results/relatedness/king/king_kinship_distribution.png",
        scatter="results/relatedness/king/king_kinship_vs_ibs0.png",
        heatmap="results/relatedness/king/king_kinship_heatmap.png"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail

        python {input.script} \
            --kin0 {input.kin0} \
            --fam {input.fam} \
            --pairs-out {output.pairs} \
            --candidates-out {output.candidates} \
            --counts-out {output.counts} \
            --summary-out {output.summary} \
            --individual-summary-out {output.individual_summary} \
            --hist-out {output.hist} \
            --scatter-out {output.scatter} \
            --heatmap-out {output.heatmap}
        """


rule write_king_relatedness_readme:
    input:
        summary="results/relatedness/king/king_relatedness_summary.tsv"
    output:
        "results/relatedness/king/README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
KING relatedness analysis
=========================

Purpose
-------
Estimate pairwise relatedness among the 46 retained Cinque Terre individuals
using the native KING-Robust kinship estimator.

KING is run with:
    king -b <PLINK BED> --kinship --cpus 1

The workflow constrains KING to one CPU for conservative laptop execution.
This computational setting does not change the estimator.

No MAF filter and no LD pruning are applied specifically for KING.

Interpretive thresholds
-----------------------
Kinship > 0.354       duplicate / monozygotic twin
0.177 to 0.354        first-degree
0.0884 to <0.177      second-degree
0.0442 to <0.0884     possible third-degree candidate
<0.0442               unrelated or more distant

For --kinship, KING documents strongest relationship-inference accuracy through
second-degree relatives. The 0.0442 third-degree boundary is retained only as
an exploratory screening threshold, not as a definitive classification.

No individual is removed automatically.

Outputs include the complete pair table, a table containing pairs crossing the
exploratory third-degree threshold, relationship-class counts, a per-individual
summary of pairwise kinship values, a kinship histogram, the recommended
kinship-versus-IBS0 diagnostic, and a pairwise kinship heatmap.

Negative KING-Robust estimates are retained as estimated. They do not represent
"negative biological relatedness" and must not be truncated to zero. Population
structure, ancestry differences and inbreeding can affect the estimator; these
features are therefore interpreted jointly with later PCA and ROH analyses.
See docs/methodological_references.md for the supporting references.

The heatmap omits the diagonal rather than assigning an artificial self-kinship
value. Its scale is therefore determined only by observed between-individual
pairwise estimates and is centered on zero when both negative and positive
estimates are present.

If close relatives are detected, sample handling is reviewed before downstream
analyses. A relatedness-driven sample change would require recalculation of
cohort-dependent site QC and HWE on the final retained sample set.
EOF
        """
