# ============================================================
# Core target QC filters before HWE
# ============================================================


rule apply_core_variant_and_individual_missingness_filters:
    input:
        vcf="results/qc/00_post_sample_qc_input/cinque_terre.post_sample_qc.autosomes.vcf.gz",
        index="results/qc/00_post_sample_qc_input/cinque_terre.post_sample_qc.autosomes.vcf.gz.tbi",
        count_script="workflow/scripts/summarize_core_filter_steps.py"
    output:
        vcf=temp("results/qc/03_core_filtering/cinque_terre.core_qc.pre_hwe.vcf.gz"),
        index=temp("results/qc/03_core_filtering/cinque_terre.core_qc.pre_hwe.vcf.gz.tbi"),
        step_counts="results/qc/03_core_filtering/core_variant_filter_step_counts.tsv",
        individual_missingness="results/qc/03_core_filtering/individual_missingness_after_site_filters.tsv",
        excluded_missingness="results/qc/03_core_filtering/individuals_excluded_missingness_gt_0.05.txt",
        retained_samples="results/qc/03_core_filtering/retained_samples_before_hwe.txt",
        summary="results/qc/03_core_filtering/core_qc_pre_hwe_summary.tsv",
        readme="results/qc/03_core_filtering/README.txt"
    log:
        "logs/qc/03_core_filtering/core_filters.log"
    params:
        variant_missingness=lambda wildcards: config["qc"]["missingness"]["remove_variants_with_missingness_gt"],
        individual_missingness=lambda wildcards: config["qc"]["missingness"]["remove_individuals_with_missingness_gt"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/03_core_filtering logs/qc/03_core_filtering

        # Exact sequential counts before filtering. The order is:
        # monomorphic -> palindromic -> variant missingness.
        {{
            printf 'REF\tALT\tMAF\tF_MISSING\n'
            bcftools +fill-tags {input.vcf} -Ou -- -t MAF,F_MISSING 2> {log}             | bcftools query -f '%REF\t%ALT\t%INFO/MAF\t%INFO/F_MISSING\n' 2>> {log}
        }} | python {input.count_script} --out {output.step_counts}

        site_tmp="results/qc/03_core_filtering/.site_filtered.tmp.vcf.gz"

        # Apply the agreed site filters. No QUAL/DP/GQ hard filters are added.
        bcftools +fill-tags {input.vcf} -Ou -- -t MAF,F_MISSING 2>> {log}         | bcftools view             -e 'MAF=0 || ((REF="A" && ALT="T") || (REF="T" && ALT="A") || (REF="C" && ALT="G") || (REF="G" && ALT="C")) || F_MISSING>{params.variant_missingness}'             -Oz -o "$site_tmp" 2>> {log}

        tabix -f -p vcf "$site_tmp" 2>> {log}
        n_site_filtered=$(bcftools index -n "$site_tmp")

        # Exact per-individual missingness after site-level filters.
        bcftools stats -s - "$site_tmp" 2>> {log}         | awk -v n_sites="$n_site_filtered" -v cutoff="{params.individual_missingness}" -F '\t' '
            BEGIN {{
                OFS="\t";
                print "sample","n_sites","n_missing","missing_rate","above_cutoff"
            }}
            $1=="PSC" {{
                rate=$14/n_sites;
                flag=(rate>cutoff ? "YES" : "NO");
                printf "%s\t%d\t%d\t%.8f\t%s\n", $3,n_sites,$14,rate,flag
            }}
        ' > {output.individual_missingness}

        awk -F '\t' 'NR>1 && $5=="YES" {{print $1}}' {output.individual_missingness}             > {output.excluded_missingness}

        n_excluded=$(wc -l < {output.excluded_missingness})

        if [ "$n_excluded" -gt 0 ]; then
            # If samples are removed, recalculate MAF and site missingness because
            # sample removal can create new monomorphic/high-missingness sites.
            bcftools view -S ^{output.excluded_missingness} "$site_tmp" -Ou 2>> {log}             | bcftools +fill-tags -Ou -- -t MAF,F_MISSING 2>> {log}             | bcftools view                 -e 'MAF=0 || F_MISSING>{params.variant_missingness}'                 -Oz -o {output.vcf} 2>> {log}
            tabix -f -p vcf {output.vcf} 2>> {log}
            rm -f "$site_tmp" "$site_tmp.tbi"
        else
            mv "$site_tmp" {output.vcf}
            mv "$site_tmp.tbi" {output.index}
        fi

        if [ ! -f {output.index} ]; then
            tabix -f -p vcf {output.vcf} 2>> {log}
        fi

        bcftools query -l {output.vcf} > {output.retained_samples}
        n_final_samples=$(wc -l < {output.retained_samples})
        n_final_variants=$(bcftools index -n {output.vcf})

        {{
            printf 'metric\tvalue\n'
            printf 'variant_missingness_cutoff_gt\t{params.variant_missingness}\n'
            printf 'individual_missingness_cutoff_gt\t{params.individual_missingness}\n'
            printf 'individuals_removed_for_missingness\t%s\n' "$n_excluded"
            printf 'retained_samples_before_hwe\t%s\n' "$n_final_samples"
            printf 'variants_before_hwe\t%s\n' "$n_final_variants"
            printf 'additional_qual_dp_gq_hard_filters\tNO\n'
        }} > {output.summary}

        cat > {output.readme} <<'EOF'
03_core_filtering
=================

Purpose
-------
Apply the agreed core QC filters after the reviewed sample exclusions.

Order of filters
----------------
1. Remove variants that became monomorphic after sample QC.
2. Remove palindromic A/T and C/G SNPs.
3. Remove variants with >5% missing genotypes.
4. Recalculate exact individual missingness and remove individuals with >5%.
5. If any individual is removed, recalculate MAF and variant missingness and
   remove newly monomorphic or >5%-missing sites.

No additional QUAL, DP, GQ or allele-balance hard filter is applied.

core_variant_filter_step_counts.tsv reports sequential counts, so overlapping
categories are not double-counted.

The VCF produced here is an intermediate pre-HWE dataset. HWE is applied in
the next rule.
EOF
        """


rule apply_hwe_bonferroni_by_chromosome:
    input:
        vcf="results/qc/03_core_filtering/cinque_terre.core_qc.pre_hwe.vcf.gz",
        index="results/qc/03_core_filtering/cinque_terre.core_qc.pre_hwe.vcf.gz.tbi",
        script="workflow/scripts/evaluate_hwe_bonferroni.py"
    output:
        vcf="results/qc/04_hwe/cinque_terre.qc_filtered.vcf.gz",
        index="results/qc/04_hwe/cinque_terre.qc_filtered.vcf.gz.tbi",
        thresholds="results/qc/04_hwe/hwe_bonferroni_thresholds.tsv",
        failed="results/qc/04_hwe/hwe_failed_variants.tsv",
        by_chromosome="results/qc/04_hwe/hwe_bonferroni_by_chromosome.tsv",
        summary="results/qc/04_hwe/final_qc_summary.tsv",
        readme="results/qc/04_hwe/README.txt"
    log:
        "logs/qc/04_hwe/hwe.log"
    params:
        alpha=lambda wildcards: config["qc"]["hwe"]["alpha"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/04_hwe logs/qc/04_hwe

        printf 'chromosome\tn_tested\talpha\tbonferroni_threshold\n' > {output.thresholds}

        bcftools index -s {input.vcf}         | sort -k1,1n         | while IFS=$'\t' read -r chrom length n_tested; do
            threshold=$(awk -v alpha="{params.alpha}" -v n="$n_tested" 'BEGIN {{printf "%.17g", alpha/n}}')
            printf '%s\t%s\t%s\t%s\n' "$chrom" "$n_tested" "{params.alpha}" "$threshold"
        done >> {output.thresholds}

        # Evaluate HWE exactly on the current retained cohort and record all
        # variants that fail the chromosome-specific Bonferroni threshold.
        {{
            printf 'CHROM\tPOS\tREF\tALT\tHWE\n'
            bcftools +fill-tags {input.vcf} -Ou -- -t HWE 2> {log}             | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/HWE\n' 2>> {log}
        }} | python {input.script}             --thresholds {output.thresholds}             --failed-out {output.failed}             --summary-out {output.by_chromosome}

        # Build one exclusion expression using the 22 chromosome-specific
        # thresholds, then apply it in one streaming pass.
        expr=""
        while IFS=$'\t' read -r chrom n_tested alpha threshold; do
            [ "$chrom" = "chromosome" ] && continue
            cond="(CHROM=\"$chrom\" && INFO/HWE<$threshold)"
            if [ -z "$expr" ]; then
                expr="$cond"
            else
                expr="$expr || $cond"
            fi
        done < {output.thresholds}

        bcftools +fill-tags {input.vcf} -Ou -- -t HWE 2>> {log}         | bcftools view -e "$expr" -Oz -o {output.vcf} 2>> {log}

        tabix -f -p vcf {output.vcf} 2>> {log}

        n_before=$(bcftools index -n {input.vcf})
        n_after=$(bcftools index -n {output.vcf})
        n_removed=$((n_before - n_after))
        n_samples=$(bcftools query -l {output.vcf} | wc -l)

        {{
            printf 'metric\tvalue\n'
            printf 'hwe_alpha\t{params.alpha}\n'
            printf 'hwe_correction\tbonferroni_by_chromosome\n'
            printf 'variants_before_hwe\t%s\n' "$n_before"
            printf 'variants_removed_hwe\t%s\n' "$n_removed"
            printf 'final_variants\t%s\n' "$n_after"
            printf 'final_samples\t%s\n' "$n_samples"
        }} > {output.summary}

        cat > {output.readme} <<'EOF'
04_hwe
======

Purpose
-------
Apply the final HWE QC rule to the post-core-QC target dataset.

HWE rule
--------
HWE is recalculated from the retained genotypes using bcftools +fill-tags.

alpha = 0.05
multiple-testing correction = Bonferroni
correction scope = separately within each chromosome

For chromosome c:
    threshold_c = 0.05 / number of variants tested on chromosome c

A variant is removed when:
    HWE p-value < threshold_c

hwe_bonferroni_by_chromosome.tsv records the number tested, threshold,
number removed and number retained for each chromosome.

hwe_failed_variants.tsv records the exact variants removed by the HWE rule.

cinque_terre.qc_filtered.vcf.gz is the final target-only QC dataset before
downstream relatedness/population-structure/selection-specific processing.
EOF
        """
