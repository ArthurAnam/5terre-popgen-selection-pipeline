# ============================================================
# Selection branch: exhaustive post-SHAPEIT2 production audit
# ============================================================

SHAPEIT2_AUDIT_DIR = SHAPEIT2_DIR + "/audit"

rule audit_selection_ct_shapeit2_phasing:
    input:
        vcfs=expand(SHAPEIT2_DIR + "/input/chr{chrom}.ct.maf005.vcf.gz", chrom=SHAPEIT2_CHROMS),
        haps=expand(SHAPEIT2_DIR + "/phased/chr{chrom}.ct.maf005.phased.haps.gz", chrom=SHAPEIT2_CHROMS),
        samples=expand(SHAPEIT2_DIR + "/phased/chr{chrom}.ct.maf005.phased.sample", chrom=SHAPEIT2_CHROMS),
        logs=expand(SHAPEIT2_DIR + "/phased/chr{chrom}.phase.log", chrom=SHAPEIT2_CHROMS),
        check=SHAPEIT2_DIR + "/check/shapeit2_check_by_chromosome.tsv",
        benchmark_haps=SHAPEIT2_DIR + "/benchmark/chr20.states400.thread1.phased.haps.gz",
        benchmark_sample=SHAPEIT2_DIR + "/benchmark/chr20.states400.thread1.phased.sample",
        script="workflow/scripts/audit_shapeit2_phasing.py"
    output:
        summary=SHAPEIT2_AUDIT_DIR + "/phasing_audit_summary.tsv",
        chromosomes=SHAPEIT2_AUDIT_DIR + "/phasing_audit_by_chromosome.tsv",
        reproducibility=SHAPEIT2_AUDIT_DIR + "/chr20_repeatability.tsv",
        sample_order=SHAPEIT2_AUDIT_DIR + "/sample_order.txt",
        gzip_ok=SHAPEIT2_AUDIT_DIR + "/gzip_integrity.ok"
    log:
        "logs/selection/phasing/shapeit2/post_phasing_audit.log"
    params:
        n=lambda wc: config["datasets"]["target"]["historical_post_qc_sample_count_to_verify"],
        ref=lambda wc: 2 * config["selection"]["phasing"]["reference_support"]["reference_samples"],
        total=lambda wc: config["selection"]["phasing"]["preflight"]["retained_snps_after_shapeit2_check"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SHAPEIT2_AUDIT_DIR} logs/selection/phasing/shapeit2
        : > {log}
        for f in {input.haps} {input.benchmark_haps}; do
            gzip -t "$f" 2>> {log}
        done
        printf 'PASS\n' > {output.gzip_ok}
        python {input.script} \
            --vcfs {input.vcfs} \
            --haps {input.haps} \
            --samples {input.samples} \
            --logs {input.logs} \
            --check {input.check} \
            --benchmark-haps {input.benchmark_haps} \
            --benchmark-sample {input.benchmark_sample} \
            --n {params.n} --ref {params.ref} --total {params.total} \
            --summary {output.summary} \
            --by-chr {output.chromosomes} \
            --repro {output.reproducibility} \
            --sample-order {output.sample_order} \
            >> {log} 2>&1
        """
