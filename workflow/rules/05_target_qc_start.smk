# ============================================================
# Start of target QC: reviewed sample removal + QC exploration
# ============================================================

QC_AUTOSOMES = ",".join(str(chrom) for chrom in range(1, 23))


rule qc_reviewed_samples_autosomes:
    input:
        vcf=lambda wildcards: config["local_paths"]["target_vcf"],
        decisions="config/sample_qc_decisions.tsv"
    output:
        excluded="results/qc/00_start/excluded_samples.txt",
        vcf="results/qc/00_start/cinque_terre.46samples.autosomes.vcf.gz",
        index="results/qc/00_start/cinque_terre.46samples.autosomes.vcf.gz.tbi",
        samples="results/qc/00_start/samples_46.txt",
        summary="results/qc/00_start/start_summary.tsv"
    log:
        "logs/qc/00_start/reviewed_samples_autosomes.log"
    params:
        regions=QC_AUTOSOMES,
        expected_samples=46
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/00_start logs/qc/00_start

        awk -F '\t' 'NR > 1 && $2 == "exclude" {{print $1}}' {input.decisions} \
            > {output.excluded}

        n_excluded=$(wc -l < {output.excluded})
        if [ "$n_excluded" -ne 4 ]; then
            echo "ERROR: expected 4 reviewed exclusions, found $n_excluded" >&2
            exit 1
        fi

        bcftools view \
            -S ^{output.excluded} \
            -r {params.regions} \
            -Oz \
            -o {output.vcf} \
            {input.vcf} 2> {log}

        tabix -f -p vcf {output.vcf} 2>> {log}
        bcftools query -l {output.vcf} > {output.samples} 2>> {log}

        n_samples=$(wc -l < {output.samples})
        n_variants=$(bcftools view -H {output.vcf} 2>> {log} | wc -l)

        if [ "$n_samples" -ne {params.expected_samples} ]; then
            echo "ERROR: expected {params.expected_samples} retained samples, found $n_samples" >&2
            exit 1
        fi

        while read -r sample; do
            if grep -Fxq "$sample" {output.samples}; then
                echo "ERROR: reviewed excluded sample still present: $sample" >&2
                exit 1
            fi
        done < {output.excluded}

        {{
            printf 'metric\tvalue\n'
            printf 'input_samples\t50\n'
            printf 'reviewed_excluded_samples\t%s\n' "$n_excluded"
            printf 'retained_samples\t%s\n' "$n_samples"
            printf 'scope\tautosomes_1_22\n'
            printf 'variants_before_variant_qc\t%s\n' "$n_variants"
            printf 'variant_filtering_applied\tNO\n'
        }} > {output.summary}
        """


rule qc_explore_site_annotations:
    input:
        vcf="results/qc/00_start/cinque_terre.46samples.autosomes.vcf.gz",
        index="results/qc/00_start/cinque_terre.46samples.autosomes.vcf.gz.tbi"
    output:
        summary="results/qc/01_exploration/site_annotation_summary.tsv",
        inventory="results/qc/01_exploration/site_inventory.tsv",
        culprit="results/qc/01_exploration/vqsr_culprit_counts.tsv",
        plot="results/qc/01_exploration/site_annotation_distributions.png"
    log:
        "logs/qc/01_exploration/site_annotations.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/01_exploration logs/qc/01_exploration

        {{
            printf 'CHROM\tPOS\tREF\tALT\tQUAL\tVQSLOD\tQD\tFS\tSOR\tMQ\tBaseQRankSum\tMQRankSum\tReadPosRankSum\tINFO_DP\tExcessHet\tInbreedingCoeff\tAC\tAN\tMAF\tF_MISSING\tculprit\tPOSITIVE_TRAIN_SITE\tNEGATIVE_TRAIN_SITE\n'

            bcftools +fill-tags {input.vcf} -Ou -- -t AC,AN,AF,MAF,F_MISSING 2> {log} \
            | bcftools query \
                -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/VQSLOD\t%INFO/QD\t%INFO/FS\t%INFO/SOR\t%INFO/MQ\t%INFO/BaseQRankSum\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\t%INFO/DP\t%INFO/ExcessHet\t%INFO/InbreedingCoeff\t%INFO/AC\t%INFO/AN\t%INFO/MAF\t%INFO/F_MISSING\t%INFO/culprit\t%INFO/POSITIVE_TRAIN_SITE\t%INFO/NEGATIVE_TRAIN_SITE\n' \
                2>> {log}
        }} \
        | python workflow/scripts/explore_variant_annotations.py \
            --summary-out {output.summary} \
            --inventory-out {output.inventory} \
            --culprit-out {output.culprit} \
            --plot-out {output.plot}
        """
