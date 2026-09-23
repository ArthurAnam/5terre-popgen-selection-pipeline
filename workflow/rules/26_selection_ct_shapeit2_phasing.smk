# ============================================================
# Selection branch: formal SHAPEIT2 alignment check and CT phasing
# ============================================================

SHAPEIT2_DIR = "results/selection/phasing/shapeit2"
SHAPEIT2_CHROMS = [str(chrom) for chrom in range(1, 23)]


def shapeit2_binary(wildcards):
    return config["local_paths"]["shapeit2_binary"]


def shapeit2_ref_haps(wildcards):
    return config["local_paths"]["shapeit2_eur_reference_haps_template"].format(chrom=wildcards.chrom)


def shapeit2_ref_legend(wildcards):
    return config["local_paths"]["shapeit2_eur_reference_legend_template"].format(chrom=wildcards.chrom)


def shapeit2_ref_sample(wildcards):
    return config["local_paths"]["shapeit2_eur_reference_sample"]


def shapeit2_map(wildcards):
    return config["local_paths"]["shapeit2_b37_map_template"].format(chrom=wildcards.chrom)


rule prepare_selection_ct_shapeit2_input:
    input:
        vcf="results/qc/04_hwe/cinque_terre.qc_filtered.vcf.gz",
        index="results/qc/04_hwe/cinque_terre.qc_filtered.vcf.gz.tbi"
    output:
        vcf=SHAPEIT2_DIR + "/input/chr{chrom}.ct.maf005.vcf.gz",
        index=SHAPEIT2_DIR + "/input/chr{chrom}.ct.maf005.vcf.gz.tbi"
    log:
        "logs/selection/phasing/shapeit2/chr{chrom}.input.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SHAPEIT2_DIR}/input logs/selection/phasing/shapeit2
        bcftools view \
            -r {wildcards.chrom} \
            -m2 -M2 -v snps \
            -Ou {input.vcf} 2> {log} \
        | bcftools +fill-tags -Ou -- -t AF 2>> {log} \
        | bcftools view \
            -i 'INFO/AF>=0.05 && INFO/AF<=0.95' \
            -Oz -o {output.vcf} 2>> {log}
        bcftools index -f -t {output.vcf} 2>> {log}
        """


rule selection_ct_shapeit2_check:
    input:
        shapeit=shapeit2_binary,
        vcf=SHAPEIT2_DIR + "/input/chr{chrom}.ct.maf005.vcf.gz",
        vcf_index=SHAPEIT2_DIR + "/input/chr{chrom}.ct.maf005.vcf.gz.tbi",
        ref_haps=shapeit2_ref_haps,
        ref_legend=shapeit2_ref_legend,
        ref_sample=shapeit2_ref_sample,
        map=shapeit2_map
    output:
        log=SHAPEIT2_DIR + "/check/chr{chrom}.check.log",
        strand=SHAPEIT2_DIR + "/check/chr{chrom}.check.snp.strand",
        exclude=SHAPEIT2_DIR + "/check/chr{chrom}.check.snp.strand.exclude"
    params:
        prefix=lambda wc: f"{SHAPEIT2_DIR}/check/chr{wc.chrom}.check"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SHAPEIT2_DIR}/check
        "{input.shapeit}" -check \
            --input-vcf {input.vcf} \
            --input-map {input.map} \
            --input-ref {input.ref_haps} {input.ref_legend} {input.ref_sample} \
            --output-log {params.prefix}

        # SHAPEIT2 normally creates both files. Keep empty files explicit if
        # a chromosome has no alignment problems so downstream rules are stable.
        [ -f {output.strand} ] || : > {output.strand}
        [ -f {output.exclude} ] || : > {output.exclude}
        """


rule summarize_selection_ct_shapeit2_check:
    input:
        preflight="results/selection/phasing/preflight/reference_overlap_by_chromosome.tsv",
        exclusions=expand(
            SHAPEIT2_DIR + "/check/chr{chrom}.check.snp.strand.exclude",
            chrom=SHAPEIT2_CHROMS,
        ),
        script="workflow/scripts/summarize_shapeit2_checks.py"
    output:
        summary=SHAPEIT2_DIR + "/check/shapeit2_check_summary.tsv",
        chromosomes=SHAPEIT2_DIR + "/check/shapeit2_check_by_chromosome.tsv"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --preflight {input.preflight} \
            --exclusions {input.exclusions} \
            --out {output.summary} \
            --chrom-out {output.chromosomes}
        """


rule audit_selection_ct_shapeit2_check:
    input:
        SHAPEIT2_DIR + "/check/shapeit2_check_summary.tsv",
        SHAPEIT2_DIR + "/check/shapeit2_check_by_chromosome.tsv"


rule selection_ct_shapeit2_phase:
    input:
        shapeit=shapeit2_binary,
        vcf=SHAPEIT2_DIR + "/input/chr{chrom}.ct.maf005.vcf.gz",
        ref_haps=shapeit2_ref_haps,
        ref_legend=shapeit2_ref_legend,
        ref_sample=shapeit2_ref_sample,
        map=shapeit2_map,
        exclude=SHAPEIT2_DIR + "/check/chr{chrom}.check.snp.strand.exclude"
    output:
        haps=SHAPEIT2_DIR + "/phased/chr{chrom}.ct.maf005.phased.haps.gz",
        sample=SHAPEIT2_DIR + "/phased/chr{chrom}.ct.maf005.phased.sample",
        log=SHAPEIT2_DIR + "/phased/chr{chrom}.phase.log"
    params:
        ne=lambda wc: config["selection"]["phasing"]["parameters"]["effective_population_size"],
        window=lambda wc: config["selection"]["phasing"]["parameters"]["window_mb"],
        states=lambda wc: config["selection"]["phasing"]["parameters"]["states"],
        burn=lambda wc: config["selection"]["phasing"]["parameters"]["burn"],
        prune=lambda wc: config["selection"]["phasing"]["parameters"]["prune"],
        main=lambda wc: config["selection"]["phasing"]["parameters"]["main"],
        seed=lambda wc: config["selection"]["phasing"]["parameters"]["seed"]
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SHAPEIT2_DIR}/phased
        "{input.shapeit}" \
            --input-vcf {input.vcf} \
            --input-map {input.map} \
            --input-ref {input.ref_haps} {input.ref_legend} {input.ref_sample} \
            --exclude-snp {input.exclude} \
            --effective-size {params.ne} \
            --window {params.window} \
            --states {params.states} \
            --burn {params.burn} \
            --prune {params.prune} \
            --main {params.main} \
            --thread {threads} \
            --seed {params.seed} \
            --output-max {output.haps} {output.sample} \
            --output-log {output.log}
        """


rule selection_ct_shapeit2_phasing:
    input:
        expand(SHAPEIT2_DIR + "/phased/chr{chrom}.ct.maf005.phased.haps.gz", chrom=SHAPEIT2_CHROMS),
        expand(SHAPEIT2_DIR + "/phased/chr{chrom}.ct.maf005.phased.sample", chrom=SHAPEIT2_CHROMS),
        expand(SHAPEIT2_DIR + "/phased/chr{chrom}.phase.log", chrom=SHAPEIT2_CHROMS)
