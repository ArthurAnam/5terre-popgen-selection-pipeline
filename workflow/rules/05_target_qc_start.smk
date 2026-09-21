# ============================================================
# Target QC start: approved sample exclusions + descriptive exploration
# ============================================================

QC_AUTOSOMES = ",".join(str(chrom) for chrom in range(1, 23))


rule build_post_sample_qc_autosomal_vcf:
    input:
        vcf=lambda wildcards: config["local_paths"]["target_vcf"],
        decisions="config/sample_qc_decisions.tsv"
    output:
        input_samples="results/qc/00_post_sample_qc_input/input_samples_before_sample_qc.txt",
        excluded="results/qc/00_post_sample_qc_input/reviewed_excluded_samples.txt",
        retained="results/qc/00_post_sample_qc_input/retained_samples_after_sample_qc.txt",
        vcf="results/qc/00_post_sample_qc_input/cinque_terre.post_sample_qc.autosomes.vcf.gz",
        index="results/qc/00_post_sample_qc_input/cinque_terre.post_sample_qc.autosomes.vcf.gz.tbi",
        summary="results/qc/00_post_sample_qc_input/post_sample_qc_input_summary.tsv",
        readme="results/qc/00_post_sample_qc_input/README.txt"
    log:
        "logs/qc/00_post_sample_qc_input/build_post_sample_qc_autosomal_vcf.log"
    params:
        regions=QC_AUTOSOMES
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/00_post_sample_qc_input logs/qc/00_post_sample_qc_input

        bcftools query -l {input.vcf} > {output.input_samples}

        awk -F '\t' 'NR > 1 && $2 == "exclude" {{print $1}}' {input.decisions} \
            > {output.excluded}

        while read -r sample; do
            if ! grep -Fxq "$sample" {output.input_samples}; then
                echo "ERROR: reviewed exclusion is not present in input VCF: $sample" >&2
                exit 1
            fi
        done < {output.excluded}

        n_input=$(wc -l < {output.input_samples})
        n_excluded=$(wc -l < {output.excluded})
        expected_retained=$((n_input - n_excluded))

        bcftools view \
            -S ^{output.excluded} \
            -r {params.regions} \
            -Oz \
            -o {output.vcf} \
            {input.vcf} 2> {log}

        tabix -f -p vcf {output.vcf} 2>> {log}
        bcftools query -l {output.vcf} > {output.retained} 2>> {log}

        n_retained=$(wc -l < {output.retained})
        n_variants=$(bcftools view -H {output.vcf} 2>> {log} | wc -l)

        if [ "$n_retained" -ne "$expected_retained" ]; then
            echo "ERROR: expected $expected_retained retained samples, found $n_retained" >&2
            exit 1
        fi

        while read -r sample; do
            if grep -Fxq "$sample" {output.retained}; then
                echo "ERROR: reviewed excluded sample still present: $sample" >&2
                exit 1
            fi
        done < {output.excluded}

        {{
            printf 'metric\tvalue\n'
            printf 'input_samples\t%s\n' "$n_input"
            printf 'reviewed_excluded_samples\t%s\n' "$n_excluded"
            printf 'retained_samples\t%s\n' "$n_retained"
            printf 'scope\tautosomes_1_22\n'
            printf 'variants_before_variant_qc\t%s\n' "$n_variants"
            printf 'variant_filtering_applied\tNO\n'
        }} > {output.summary}

        cat > {output.readme} <<'EOF'
00_post_sample_qc_input
=======================

Purpose
-------
This folder contains the official starting dataset for variant-level QC.

No variant-level quality filter has been applied here.

Files
-----
input_samples_before_sample_qc.txt
    All sample IDs present in the delivered input VCF before the approved
    sample-level exclusions.

reviewed_excluded_samples.txt
    Sample IDs excluded after explicit review of the sample-level QC evidence.
    "Reviewed" means the exclusion was manually/scientifically approved and was
    not triggered automatically by a single threshold.

retained_samples_after_sample_qc.txt
    Sample IDs retained after applying the reviewed exclusions.

cinque_terre.post_sample_qc.autosomes.vcf.gz
    Cinque Terre VCF after the approved sample exclusions and restriction to
    autosomes 1-22. This is the starting VCF for variant-level QC.

cinque_terre.post_sample_qc.autosomes.vcf.gz.tbi
    Tabix index for the VCF above.

post_sample_qc_input_summary.tsv
    Counts of input samples, excluded samples, retained samples and variants.

Important
---------
This step does NOT repeat variant calling or VQSR and does NOT apply DP, GQ,
QUAL, HWE, missingness, palindromic or other variant-level filters.
EOF
        """


rule explore_variant_qc_metrics_without_filtering:
    input:
        vcf="results/qc/00_post_sample_qc_input/cinque_terre.post_sample_qc.autosomes.vcf.gz",
        index="results/qc/00_post_sample_qc_input/cinque_terre.post_sample_qc.autosomes.vcf.gz.tbi",
        script="workflow/exploratory/explore_variant_qc_metrics.py"
    output:
        summary="results/qc/01_variant_qc_exploration/variant_qc_metric_descriptive_statistics.tsv",
        inventory="results/qc/01_variant_qc_exploration/variant_category_counts_before_filtering.tsv",
        plot="results/qc/01_variant_qc_exploration/variant_qc_metric_distributions.png",
        readme="results/qc/01_variant_qc_exploration/README.txt"
    log:
        "logs/qc/01_variant_qc_exploration/variant_qc_exploration.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/01_variant_qc_exploration logs/qc/01_variant_qc_exploration

        {{
            printf 'CHROM\tPOS\tREF\tALT\tQUAL\tMAF\tF_MISSING\n'

            bcftools +fill-tags {input.vcf} -Ou -- -t AC,AN,AF,MAF,F_MISSING 2> {log} \
            | bcftools query \
                -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/MAF\t%INFO/F_MISSING\n' \
                2>> {log}
        }} \
        | python {input.script} \
            --summary-out {output.summary} \
            --inventory-out {output.inventory} \
            --plot-out {output.plot}

        cat > {output.readme} <<'EOF'
01_variant_qc_exploration
=========================

Purpose
-------
Descriptive exploration before variant filtering.

Nothing in this folder is an automatic filtering rule.

Current metrics
---------------
QUAL
    VCF site-level variant quality score carried by the delivered callset.
    It is inspected descriptively only.

MAF
    Minor allele frequency recalculated on the retained post-sample-QC cohort.

F_MISSING
    Fraction of missing genotypes recalculated on the retained cohort.

variant_category_counts_before_filtering.tsv also reports:
    total sites;
    biallelic SNPs;
    multiallelic/non-SNP records;
    variants that became monomorphic after sample QC;
    palindromic A/T or C/G SNPs;
    variants with >5% missing genotypes.

Why VQSR fields are not included here
--------------------------------------
VQSR was completed upstream before this pipeline. VQSLOD, culprit,
POSITIVE_TRAIN_SITE and NEGATIVE_TRAIN_SITE describe the upstream VQSR model
and are not part of the variant-level QC decisions performed here.

Likewise, this stage does not attempt to reproduce or second-guess VQSR.

Next diagnostic block
---------------------
FORMAT-level DP, GQ and AD will be explored separately because they describe
individual genotype calls rather than site-level VQSR bookkeeping. They will
first be inspected without hard filtering.
EOF
        """
