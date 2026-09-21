# ============================================================
# KING relatedness analysis: Cinque Terre only
# ============================================================

rule prepare_king_relatedness_dataset:
    input:
        vcf="results/qc/03_core_filtering/cinque_terre.core_qc.pre_hwe.vcf.gz",
        index="results/qc/03_core_filtering/cinque_terre.core_qc.pre_hwe.vcf.gz.tbi"
    output:
        bed=temp("results/relatedness/king/input/cinque_terre.king.bed"),
        bim=temp("results/relatedness/king/input/cinque_terre.king.bim"),
        fam=temp("results/relatedness/king/input/cinque_terre.king.fam"),
        summary="results/relatedness/king/king_input_dataset_summary.tsv",
        readme="results/relatedness/king/INPUT_README.txt"
    log:
        "logs/relatedness/king/prepare_king_input.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/relatedness/king/input logs/relatedness/king

        plink2 \
            --vcf {input.vcf} \
            --double-id \
            --set-all-var-ids '@:#:$r:$a' \
            --make-bed \
            --out results/relatedness/king/input/cinque_terre.king \
            > {log} 2>&1

        n_samples=$(wc -l < {output.fam})
        n_variants=$(wc -l < {output.bim})

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
        }} > {output.summary}

        cat > {output.readme} <<'EOF'
KING input dataset
==================

The KING dataset is derived from the stable Cinque Terre QC dataset before HWE filtering.

No extra population-genetic filter is added specifically for KING:
- no MAF filter;
- no LD pruning.

HWE filtering is intentionally not used to define the KING input. Relatedness
is being evaluated as a QC/property of the retained cohort, and KING recommends
retaining genome-wide SNPs that already passed technical QC rather than pruning
or filtering otherwise good SNPs.

PLINK2 converts the stable pre-HWE QC VCF to PLINK BED/BIM/FAM format. Each individual
is assigned a unique family ID (--double-id), so all 46 choose 2 pairwise
comparisons are treated as between-family comparisons and are expected in
KING's .kin0 output.

The PLINK binary files are temporary workflow products and can be regenerated
from the stable pre-HWE QC VCF.
EOF
        """

rule run_king_relatedness:
    input:
        bed="results/relatedness/king/input/cinque_terre.king.bed",
        bim="results/relatedness/king/input/cinque_terre.king.bim",
        fam="results/relatedness/king/input/cinque_terre.king.fam",
        script="workflow/scripts/plot_king_relatedness.py"
    output:
        kin0="results/relatedness/king/king.kin0",
        pairs="results/relatedness/king/king_pairwise_kinship.tsv",
        candidates="results/relatedness/king/king_candidate_relatives.tsv",
        counts="results/relatedness/king/king_relationship_counts.tsv",
        summary="results/relatedness/king/king_relatedness_summary.tsv",
        hist="results/relatedness/king/king_kinship_distribution.png",
        scatter="results/relatedness/king/king_kinship_vs_ibs0.png",
        heatmap="results/relatedness/king/king_kinship_heatmap.png",
        readme="results/relatedness/king/README.txt"
    log:
        "logs/relatedness/king/king.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/relatedness/king logs/relatedness/king

        king \
            -b {input.bed} \
            --kinship \
            --prefix results/relatedness/king/king \
            > {log} 2>&1

        test -s {output.kin0}

        python {input.script} \
            --kin0 {output.kin0} \
            --fam {input.fam} \
            --pairs-out {output.pairs} \
            --candidates-out {output.candidates} \
            --counts-out {output.counts} \
            --summary-out {output.summary} \
            --hist-out {output.hist} \
            --scatter-out {output.scatter} \
            --heatmap-out {output.heatmap}

        cat > {output.readme} <<'EOF'
KING relatedness analysis
=========================

Purpose
-------
Estimate pairwise relatedness among the 46 retained Cinque Terre individuals
using the native KING robust kinship estimator.

KING is run with:
    king -b <PLINK BED> --kinship

No MAF filter and no LD pruning are applied specifically for KING.

Interpretive thresholds
-----------------------
Kinship > 0.354       duplicate / monozygotic twin
0.177 to 0.354        first-degree
0.0884 to <0.177      second-degree
0.0442 to <0.0884     possible third-degree candidate
<0.0442               unrelated or more distant

For --kinship, KING documents strongest relationship-inference accuracy through
second-degree relatives. The 0.0442 third-degree boundary is therefore retained
only as an exploratory screening threshold, not as a definitive classification.

No individual is removed automatically.

Outputs include the complete pair table, a table of pairs at third degree or
closer, relationship-class counts, a kinship histogram, the classic
kinship-versus-IBS0 diagnostic, and a pairwise kinship heatmap.

If close relatives are detected, sample handling will be reviewed before
downstream analyses. Since HWE has already been applied in the current QC,
a substantive relatedness finding would also motivate revisiting HWE after
the final sample set is decided.
EOF
        """
