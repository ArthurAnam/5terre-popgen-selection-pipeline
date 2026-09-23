# ============================================================
# Selection branch: CT SHAPEIT2 EUR-reference alignment preflight
# ============================================================

PHASING_PREFLIGHT_DIR = "results/selection/phasing/preflight"
PHASING_CHROMS = [str(chrom) for chrom in range(1, 23)]


rule selection_ct_phasing_overlap_audit_chromosome:
    input:
        ct_vcf="results/qc/04_hwe/cinque_terre.qc_filtered.vcf.gz",
        ct_index="results/qc/04_hwe/cinque_terre.qc_filtered.vcf.gz.tbi",
        ct_maf005_bim=SEL_LD_DIR + "/panels/maf005/ct.maf005.bim",
        eur_samples="results/reference/1kg_eur/metadata/eur_sample_ids.txt",
        ref_vcf=onekg_phase3_vcf,
        ref_index=onekg_phase3_vcf_index
    output:
        audit=PHASING_PREFLIGHT_DIR + "/per_chromosome/chr{chrom}.reference_overlap.tsv"
    log:
        "logs/selection/phasing/preflight/chr{chrom}.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        export LC_ALL=C
        mkdir -p {PHASING_PREFLIGHT_DIR}/per_chromosome logs/selection/phasing/preflight
        : > {log}

        tmpdir="{PHASING_PREFLIGHT_DIR}/.tmp_chr{wildcards.chrom}.$$"
        mkdir -p "$tmpdir"
        trap 'rm -rf "$tmpdir"' EXIT

        ct="$tmpdir/ct.maf005.bcf"
        ref_eur="$tmpdir/ref.eur.poly.bcf"
        shared_eur="$tmpdir/shared.eur.bcf"

        bcftools view \
            -r {wildcards.chrom} \
            -m2 -M2 -v snps \
            -Ou {input.ct_vcf} 2>> {log} \
        | bcftools +fill-tags -Ou -- -t AC,AN,AF 2>> {log} \
        | bcftools view \
            -i 'INFO/AF>=0.05 && INFO/AF<=0.95' \
            -Ob -o "$ct" 2>> {log}
        bcftools index -f "$ct" 2>> {log}

        bcftools view \
            -S {input.eur_samples} \
            -m2 -M2 -v snps \
            -Ou {input.ref_vcf} 2>> {log} \
        | bcftools +fill-tags -Ou -- -t AC,AN 2>> {log} \
        | bcftools view \
            -i 'INFO/AC>0 && INFO/AC<INFO/AN' \
            -Ob -o "$ref_eur" 2>> {log}
        bcftools index -f "$ref_eur" 2>> {log}

        bcftools isec -c none -n=2 -w1 \
            -Ob -o "$shared_eur" "$ct" "$ref_eur" 2>> {log}
        bcftools index -f "$shared_eur" 2>> {log}

        bcftools query -f '%CHROM\t%POS\n' "$ct" | sort -u > "$tmpdir/ct.pos"
        bcftools query -f '%CHROM\t%POS\n' "$ref_eur" | sort -u > "$tmpdir/eur.pos"

        ct_n=$(bcftools index -n "$ct")
        ct_plink=$(awk -v c="{wildcards.chrom}" '$1==c {{n++}} END {{print n+0}}' {input.ct_maf005_bim})
        eur_n=$(bcftools index -n "$ref_eur")
        eur_exact=$(bcftools index -n "$shared_eur")
        eur_same=$(comm -12 "$tmpdir/ct.pos" "$tmpdir/eur.pos" | wc -l)
        eur_samples=$(bcftools query -l "$ref_eur" | wc -l)

        if [ "$ct_n" -eq "$ct_plink" ]; then panel_match="PASS"; else panel_match="FAIL"; fi

        printf 'chromosome\tct_maf005_bcftools\tct_maf005_plink\tct_panel_count_match\teur_polymorphic_biallelic_snps\teur_samples\texact_eur_polymorphic\tfraction_ct_exact_eur_polymorphic\tsame_position_eur_polymorphic\tallele_mismatch_eur_polymorphic\n' \
            > {output.audit}
        awk -v chr="{wildcards.chrom}" \
            -v ct="$ct_n" -v ctp="$ct_plink" -v pm="$panel_match" \
            -v en="$eur_n" -v es="$eur_samples" -v ee="$eur_exact" -v esp="$eur_same" \
            'BEGIN {{ OFS="\t"; print chr,ct,ctp,pm,en,es,ee,(ct?ee/ct:0),esp,(esp-ee) }}' \
            >> {output.audit}
        """


rule summarize_selection_ct_phasing_preflight:
    input:
        audits=expand(
            PHASING_PREFLIGHT_DIR + "/per_chromosome/chr{chrom}.reference_overlap.tsv",
            chrom=PHASING_CHROMS,
        ),
        script="workflow/scripts/summarize_phasing_preflight.py"
    output:
        summary=PHASING_PREFLIGHT_DIR + "/reference_overlap_summary.tsv",
        chromosomes=PHASING_PREFLIGHT_DIR + "/reference_overlap_by_chromosome.tsv"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --audits {input.audits} \
            --summary-out {output.summary} \
            --chrom-out {output.chromosomes}
        """


rule audit_selection_ct_phasing_preflight:
    input:
        PHASING_PREFLIGHT_DIR + "/reference_overlap_summary.tsv",
        PHASING_PREFLIGHT_DIR + "/reference_overlap_by_chromosome.tsv"
