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
        iteration_history="results/qc/03_core_filtering/iterative_missingness_qc_history.tsv",
        individual_missingness="results/qc/03_core_filtering/individual_missingness_after_site_filters.tsv",
        excluded_missingness="results/qc/03_core_filtering/individuals_excluded_missingness_gt_0.05.txt",
        retained_samples="results/qc/03_core_filtering/retained_samples_before_hwe.txt",
        summary="results/qc/03_core_filtering/core_qc_pre_hwe_summary.tsv"
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
        : > {log}

        # Audit the first sequential pass without changing the input:
        # static palindromic filter, then dynamic monomorphism and site missingness.
        {{
            printf 'REF\tALT\tAC\tAN\tF_MISSING\n'
            bcftools +fill-tags {input.vcf} -Ou -- -t AC,AN,F_MISSING 2>> {log} \
            | bcftools query -f '%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/F_MISSING\n' 2>> {log}
        }} | python {input.count_script} --out {output.step_counts}

        # Static filter: REF/ALT class does not depend on cohort composition.
        static_tmp="results/qc/03_core_filtering/.static_no_palindromic.vcf.gz"
        bcftools view \
            -e '((REF="A" && ALT="T") || (REF="T" && ALT="A") || (REF="C" && ALT="G") || (REF="G" && ALT="C"))' \
            {input.vcf} -Oz -o "$static_tmp" 2>> {log}
        tabix -f -p vcf "$static_tmp" 2>> {log}

        # Dynamic QC is monotonic and is repeated only when an individual is
        # removed. Site metrics are then recomputed on the updated cohort.
        printf 'iteration\tsamples_before\tvariants_before\tvariants_after_site_filters\tsites_removed_this_iteration\tindividuals_removed_this_iteration\tsamples_after\tconverged\n' \
            > {output.iteration_history}
        : > {output.excluded_missingness}

        current="$static_tmp"
        iteration=1
        max_iterations=50

        while true; do
            if [ "$iteration" -gt "$max_iterations" ]; then
                echo "ERROR: iterative QC did not converge within $max_iterations iterations" >&2
                exit 1
            fi

            samples_before=$(bcftools query -l "$current" | wc -l)
            variants_before=$(bcftools index -n "$current")

            site_tmp="results/qc/03_core_filtering/.dynamic_sites.iter_$iteration.vcf.gz"

            # Recompute AC/AN and site missingness on the current cohort.
            # Monomorphic means AC=0 or AC=AN with AN>0; AN=0 is removed by
            # the >5% missingness criterion rather than called monomorphic.
            bcftools +fill-tags "$current" -Ou -- -t AC,AN,F_MISSING 2>> {log} \
            | bcftools view \
                -e '((AN>0 && (AC=0 || AC=AN)) || F_MISSING>{params.variant_missingness})' \
                -Oz -o "$site_tmp" 2>> {log}
            tabix -f -p vcf "$site_tmp" 2>> {log}

            variants_after=$(bcftools index -n "$site_tmp")
            sites_removed=$((variants_before - variants_after))

            if [ "$variants_after" -eq 0 ]; then
                echo "ERROR: no variants remain after dynamic site QC" >&2
                exit 1
            fi

            miss_tmp="results/qc/03_core_filtering/.individual_missingness.iter_$iteration.tsv"
            bcftools stats -s - "$site_tmp" 2>> {log} \
            | awk -v n_sites="$variants_after" -v cutoff="{params.individual_missingness}" -F '\t' '
                BEGIN {{
                    OFS="\t";
                    print "sample","n_sites","n_missing","missing_rate","above_cutoff"
                }}
                $1=="PSC" {{
                    rate=$14/n_sites;
                    flag=(rate>cutoff ? "YES" : "NO");
                    printf "%s\t%d\t%d\t%.8f\t%s\n", $3,n_sites,$14,rate,flag
                }}
            ' > "$miss_tmp"

            excluded_iter="results/qc/03_core_filtering/.excluded.iter_$iteration.txt"
            awk -F '\t' 'NR>1 && $5=="YES" {{print $1}}' "$miss_tmp" > "$excluded_iter"
            n_excluded=$(wc -l < "$excluded_iter")
            samples_after=$((samples_before - n_excluded))

            if [ "$n_excluded" -eq 0 ]; then
                printf '%s\t%s\t%s\t%s\t%s\t0\t%s\tYES\n' \
                    "$iteration" "$samples_before" "$variants_before" "$variants_after" "$sites_removed" "$samples_after" \
                    >> {output.iteration_history}

                mv "$site_tmp" {output.vcf}
                mv "$site_tmp.tbi" {output.index}
                mv "$miss_tmp" {output.individual_missingness}
                rm -f "$excluded_iter" "$current" "$current.tbi"
                break
            fi

            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tNO\n' \
                "$iteration" "$samples_before" "$variants_before" "$variants_after" "$sites_removed" "$n_excluded" "$samples_after" \
                >> {output.iteration_history}

            cat "$excluded_iter" >> {output.excluded_missingness}
            sort -u -o {output.excluded_missingness} {output.excluded_missingness}

            if [ "$samples_after" -le 0 ]; then
                echo "ERROR: individual-missingness QC would remove all samples" >&2
                exit 1
            fi

            next_tmp="results/qc/03_core_filtering/.samples.iter_$iteration.vcf.gz"
            bcftools view -S ^"$excluded_iter" "$site_tmp" -Oz -o "$next_tmp" 2>> {log}
            tabix -f -p vcf "$next_tmp" 2>> {log}

            rm -f "$site_tmp" "$site_tmp.tbi" "$miss_tmp" "$excluded_iter" "$current" "$current.tbi"
            current="$next_tmp"
            iteration=$((iteration + 1))
        done

        bcftools query -l {output.vcf} > {output.retained_samples}
        n_final_samples=$(wc -l < {output.retained_samples})
        n_final_variants=$(bcftools index -n {output.vcf})
        n_total_excluded=$(wc -l < {output.excluded_missingness})
        n_iterations=$(awk 'END {{print NR-1}}' {output.iteration_history})

        {{
            printf 'metric\tvalue\n'
            printf 'variant_missingness_cutoff_gt\t{params.variant_missingness}\n'
            printf 'individual_missingness_cutoff_gt\t{params.individual_missingness}\n'
            printf 'dynamic_qc_iterations\t%s\n' "$n_iterations"
            printf 'individuals_removed_for_missingness\t%s\n' "$n_total_excluded"
            printf 'retained_samples_before_hwe\t%s\n' "$n_final_samples"
            printf 'variants_before_hwe\t%s\n' "$n_final_variants"
            printf 'additional_qual_dp_gq_hard_filters\tNO\n'
        }} > {output.summary}

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
        summary="results/qc/04_hwe/final_qc_summary.tsv"
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

        bcftools index -s {input.vcf} \
        | sort -k1,1n \
        | while IFS=$'\t' read -r chrom length n_tested; do
            threshold=$(awk -v alpha="{params.alpha}" -v n="$n_tested" 'BEGIN {{printf "%.17g", alpha/n}}')
            printf '%s\t%s\t%s\t%s\n' "$chrom" "$n_tested" "{params.alpha}" "$threshold"
        done >> {output.thresholds}

        # Evaluate HWE exactly on the current retained cohort and record all
        # variants that fail the chromosome-specific Bonferroni threshold.
        {{
            printf 'CHROM\tPOS\tREF\tALT\tHWE\n'
            bcftools +fill-tags {input.vcf} -Ou -- -t HWE 2> {log} \
            | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/HWE\n' 2>> {log}
        }} | python {input.script} \
            --thresholds {output.thresholds} \
            --failed-out {output.failed} \
            --summary-out {output.by_chromosome}

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

        bcftools +fill-tags {input.vcf} -Ou -- -t HWE 2>> {log} \
        | bcftools view -e "$expr" -Oz -o {output.vcf} 2>> {log}

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

        """


# Documentation-only rules are deliberately separated from the expensive
# computational rules above. Editing explanatory text must not trigger a
# full re-run of core QC or HWE.

rule write_core_qc_readme:
    input:
        summary="results/qc/03_core_filtering/core_qc_pre_hwe_summary.tsv"
    output:
        "results/qc/03_core_filtering/README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
03_core_filtering
=================

Purpose
-------
Apply the agreed core QC filters after the four reviewed sample exclusions.

Static filter
-------------
Palindromic A/T and C/G SNPs are removed first. This REF/ALT property does not
change when samples or other variants are removed.

Dynamic QC until stability
--------------------------
The following cohort-dependent quantities are then recalculated on the current
dataset:
1. monomorphic status from AC/AN;
2. variant missingness;
3. individual missingness.

Monomorphic sites and variants with >5% missing genotypes are removed, then
individual missingness is calculated on the retained sites. Individuals with
>5% missingness are removed.

If at least one individual is removed, AC/AN and variant missingness are
recalculated on the new cohort and the dynamic block is repeated. The block
stops when an iteration removes no additional individual. At that point,
repeating the site filters on the unchanged cohort cannot remove additional
sites, so the dataset is stable for these criteria.

The procedure is monotonic: once a site or individual is removed, it is not
reintroduced in a later iteration.

Audit files
-----------
core_variant_filter_step_counts.tsv reports the exact sequential effect of the
first pass: palindromic -> monomorphic -> variant missingness.

iterative_missingness_qc_history.tsv records every dynamic iteration, including
sample and variant counts and whether convergence was reached.

individual_missingness_after_site_filters.tsv contains the exact individual
missingness values from the stable iteration.

No additional QUAL, DP, GQ or allele-balance hard filter is applied.

The resulting VCF is the stable pre-HWE dataset. HWE is applied only afterward.
EOF
        """


rule write_hwe_readme:
    input:
        summary="results/qc/04_hwe/final_qc_summary.tsv"
    output:
        "results/qc/04_hwe/README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
04_hwe
======

Purpose
-------
Apply the final HWE QC rule to the stable post-core-QC target dataset.

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

cinque_terre.qc_filtered.vcf.gz is the final target-only QC dataset for
downstream population-structure and selection-specific processing.

The KING relatedness QC branch is intentionally evaluated on the stable
pre-HWE cohort so that any relatedness-driven sample decision can be made
before final HWE filtering is considered definitive.
EOF
        """
