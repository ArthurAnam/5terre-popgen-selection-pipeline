# ============================================================
# Pre-QC for the Cinque Terre target VCF
# ============================================================
# This module validates and inventories the delivered post-VQSR VCF.
# It does NOT filter or rewrite variants.

rule preqc_target_vcf:
    input:
        vcf=lambda wildcards: config["local_paths"]["target_vcf"],
        index=lambda wildcards: config["local_paths"]["target_vcf_index"]
    output:
        quickcheck="results/preqc/input.quickcheck.txt",
        samples="results/preqc/samples.txt",
        filters="results/preqc/filter_counts.tsv",
        contigs="results/preqc/contigs.txt",
        stats="results/preqc/input.bcftools.stats.txt",
        summary="results/preqc/preqc_summary.tsv",
        done="results/preqc/preqc_done.txt"
    log:
        "logs/preqc/preqc_target_vcf.log"
    params:
        expected_samples=config["datasets"]["target"]["sample_count"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/preqc logs/preqc

        # 1. Basic structural integrity check. Silence means success.
        if bcftools quickcheck -v {input.vcf} > {output.quickcheck} 2>&1; then
            :
        else
            cat {output.quickcheck} >&2
            exit 1
        fi

        # 2. The tabix index must be usable; listing indexed contigs also records naming.
        tabix -l {input.vcf} > {output.contigs} 2>> {log}

        # 3. Inventory samples and FILTER values without modifying the input VCF.
        bcftools query -l {input.vcf} > {output.samples} 2>> {log}
        bcftools query -f '%FILTER\n' {input.vcf} 2>> {log} \
            | sort \
            | uniq -c \
            | awk 'BEGIN{{OFS="\t"; print "count","FILTER"}} {{print $1,$2}}' \
            > {output.filters}

        # 4. Full baseline statistics for the delivered file.
        bcftools stats {input.vcf} > {output.stats} 2>> {log}

        # 5. Assert the two upstream facts currently expected for this dataset:
        #    46 samples and only VQSR-passing records (FILTER=PASS).
        n_samples=$(wc -l < {output.samples})
        n_variants=$(bcftools view -H {input.vcf} 2>> {log} | wc -l)
        n_nonpass=$(bcftools query -f '%FILTER\n' {input.vcf} 2>> {log} \
            | awk '$1 != "PASS" {{n++}} END {{print n+0}}')

        {{
            printf 'metric\tvalue\n'
            printf 'expected_samples\t%s\n' '{params.expected_samples}'
            printf 'observed_samples\t%s\n' "$n_samples"
            printf 'total_variants\t%s\n' "$n_variants"
            printf 'non_PASS_variants\t%s\n' "$n_nonpass"
        }} > {output.summary}

        if [ "$n_samples" -ne {params.expected_samples} ]; then
            echo "ERROR: expected {params.expected_samples} samples, found $n_samples" >&2
            exit 1
        fi

        if [ "$n_nonpass" -ne 0 ]; then
            echo "ERROR: found $n_nonpass variants whose FILTER value is not PASS" >&2
            exit 1
        fi

        printf 'Pre-QC passed: %s samples; %s variants; all FILTER=PASS.\n' \
            "$n_samples" "$n_variants" > {output.done}
        """
