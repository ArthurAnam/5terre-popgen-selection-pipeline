# ============================================================
# 1000 Genomes EUR QC + exact harmonization with Cinque Terre
# ============================================================
#
# This branch deliberately keeps the large BCF intermediates temporary.
# The persistent master output is a PLINK2 PGEN/PVAR/PSAM dataset containing
# only exact CHR:POS:REF:ALT matches between the final Cinque Terre QC dataset
# and the QCed 503-sample 1000 Genomes EUR subset.

HARMONIZATION_CHROMS = [str(chrom) for chrom in range(1, 23)]


def onekg_phase3_vcf(wildcards):
    return os.path.join(
        config["local_paths"]["onekg_phase3_vcf_dir"],
        config["datasets"]["reference_1kg_eur"]["phase3_vcf_filename_template"].format(
            chrom=wildcards.chrom
        ),
    )


def onekg_phase3_vcf_index(wildcards):
    return onekg_phase3_vcf(wildcards) + ".tbi"


rule prepare_1kg_eur_metadata:
    input:
        metadata=lambda wildcards: config["local_paths"]["onekg_phase3_sample_metadata"]
    output:
        metadata="results/reference/1kg_eur/metadata/eur_samples.tsv",
        samples="results/reference/1kg_eur/metadata/eur_sample_ids.txt",
        summary="results/reference/1kg_eur/metadata/eur_sample_inventory.tsv"
    params:
        expected_total=lambda wildcards: config["datasets"]["reference_1kg_eur"]["populations"]["expected_total_eur_samples"],
        ceu=lambda wildcards: config["datasets"]["reference_1kg_eur"]["populations"]["expected_counts"]["CEU"],
        fin=lambda wildcards: config["datasets"]["reference_1kg_eur"]["populations"]["expected_counts"]["FIN"],
        gbr=lambda wildcards: config["datasets"]["reference_1kg_eur"]["populations"]["expected_counts"]["GBR"],
        ibs=lambda wildcards: config["datasets"]["reference_1kg_eur"]["populations"]["expected_counts"]["IBS"],
        tsi=lambda wildcards: config["datasets"]["reference_1kg_eur"]["populations"]["expected_counts"]["TSI"]
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/reference/1kg_eur/metadata

        awk '
        BEGIN {{OFS="\t"}}
        NR==1 {{print "ID","POP","GROUP","SEX"; next}}
        $3=="EUR" {{print $1,$2,$3,$4}}
        ' {input.metadata} > {output.metadata}

        tail -n +2 {output.metadata} | cut -f1 > {output.samples}

        n_total=$(wc -l < {output.samples})
        if [ "$n_total" -ne "{params.expected_total}" ]; then
            echo "ERROR: expected {params.expected_total} EUR samples, found $n_total" >&2
            exit 1
        fi

        n_ceu=$(awk -F '\t' '$2=="CEU" {{n++}} END {{print n+0}}' {output.metadata})
        n_fin=$(awk -F '\t' '$2=="FIN" {{n++}} END {{print n+0}}' {output.metadata})
        n_gbr=$(awk -F '\t' '$2=="GBR" {{n++}} END {{print n+0}}' {output.metadata})
        n_ibs=$(awk -F '\t' '$2=="IBS" {{n++}} END {{print n+0}}' {output.metadata})
        n_tsi=$(awk -F '\t' '$2=="TSI" {{n++}} END {{print n+0}}' {output.metadata})

        [ "$n_ceu" -eq "{params.ceu}" ] || {{ echo "ERROR: CEU count mismatch" >&2; exit 1; }}
        [ "$n_fin" -eq "{params.fin}" ] || {{ echo "ERROR: FIN count mismatch" >&2; exit 1; }}
        [ "$n_gbr" -eq "{params.gbr}" ] || {{ echo "ERROR: GBR count mismatch" >&2; exit 1; }}
        [ "$n_ibs" -eq "{params.ibs}" ] || {{ echo "ERROR: IBS count mismatch" >&2; exit 1; }}
        [ "$n_tsi" -eq "{params.tsi}" ] || {{ echo "ERROR: TSI count mismatch" >&2; exit 1; }}

        {{
            printf 'population\tn_samples\n'
            printf 'CEU\t%s\n' "$n_ceu"
            printf 'FIN\t%s\n' "$n_fin"
            printf 'GBR\t%s\n' "$n_gbr"
            printf 'IBS\t%s\n' "$n_ibs"
            printf 'TSI\t%s\n' "$n_tsi"
            printf 'EUR_TOTAL\t%s\n' "$n_total"
        }} > {output.summary}
        """


rule harmonize_ct_1kg_eur_chromosome:
    input:
        ct_vcf="results/qc/04_hwe/cinque_terre.qc_filtered.vcf.gz",
        ct_index="results/qc/04_hwe/cinque_terre.qc_filtered.vcf.gz.tbi",
        eur_samples="results/reference/1kg_eur/metadata/eur_sample_ids.txt",
        eur_vcf=onekg_phase3_vcf,
        eur_index=onekg_phase3_vcf_index
    output:
        bcf=temp("results/population_structure/harmonization/per_chromosome/chr{chrom}.harmonized.bcf"),
        index=temp("results/population_structure/harmonization/per_chromosome/chr{chrom}.harmonized.bcf.csi"),
        audit="results/population_structure/harmonization/audit/chr{chrom}.harmonization.tsv"
    params:
        variant_missingness=lambda wildcards: config["datasets"]["reference_1kg_eur"]["downstream_qc"]["remove_variants_with_missingness_gt"]
    log:
        "logs/population_structure/harmonization/chr{chrom}.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p             results/population_structure/harmonization/per_chromosome             results/population_structure/harmonization/audit             logs/population_structure/harmonization

        : > {log}

        tmpdir="results/population_structure/harmonization/.tmp_chr{wildcards.chrom}.$$"
        mkdir -p "$tmpdir"
        trap 'rm -rf "$tmpdir"' EXIT

        ct_chr="$tmpdir/ct.bcf"
        eur_chr="$tmpdir/eur_qc.bcf"
        ct_common="$tmpdir/ct_common.bcf"
        eur_common="$tmpdir/eur_common.bcf"

        bcftools view             -r {wildcards.chrom}             -Ob -o "$ct_chr"             {input.ct_vcf} 2>> {log}
        bcftools index -f "$ct_chr" 2>> {log}

        bcftools view             -S {input.eur_samples}             -m2 -M2             -v snps             -Ou {input.eur_vcf} 2>> {log}         | bcftools +fill-tags -Ou -- -t AC,AN,F_MISSING 2>> {log}         | bcftools view             -e 'AN=0 || (AN>0 && (AC=0 || AC=AN)) || F_MISSING>{params.variant_missingness} || ((REF="A" && ALT="T") || (REF="T" && ALT="A") || (REF="C" && ALT="G") || (REF="G" && ALT="C"))'             -Ob -o "$eur_chr" 2>> {log}
        bcftools index -f "$eur_chr" 2>> {log}

        # The audited Phase 3 EUR source has zero missing genotypes after the
        # site filters above. Enforce that observation here: if a future input
        # differs, stop rather than silently skipping the planned global
        # individual-missingness QC.
        n_any_missing=$(bcftools query -f '%INFO/F_MISSING\n' "$eur_chr" 2>> {log}             | awk '$1>0 {{n++}} END {{print n+0}}')
        if [ "$n_any_missing" -ne 0 ]; then
            echo "ERROR: EUR reference contains retained sites with missing genotypes; global individual-missingness QC must be implemented before harmonization." >&2
            exit 1
        fi

        bcftools query -l "$ct_chr" | sort > "$tmpdir/ct.samples"
        bcftools query -l "$eur_chr" | sort > "$tmpdir/eur.samples"
        comm -12 "$tmpdir/ct.samples" "$tmpdir/eur.samples" > "$tmpdir/sample_overlap.txt"
        if [ -s "$tmpdir/sample_overlap.txt" ]; then
            echo "ERROR: sample IDs overlap between Cinque Terre and 1000G EUR" >&2
            cat "$tmpdir/sample_overlap.txt" >&2
            exit 1
        fi

        bcftools query -f '%CHROM\t%POS\n' "$ct_chr" 2>> {log}             | sort -k1,1 -k2,2n -u > "$tmpdir/ct.pos"
        bcftools query -f '%CHROM\t%POS\n' "$eur_chr" 2>> {log}             | sort -k1,1 -k2,2n -u > "$tmpdir/eur.pos"

        same_position=$(comm -12 "$tmpdir/ct.pos" "$tmpdir/eur.pos" | wc -l)

        # -c none requires allele-exact compatibility. Two passes retain the
        # matching records from each source separately for the subsequent
        # sample-wise merge.
        bcftools isec             -c none -n=2 -w1             -Ob -o "$ct_common"             "$ct_chr" "$eur_chr" 2>> {log}
        bcftools index -f "$ct_common" 2>> {log}

        bcftools isec             -c none -n=2 -w2             -Ob -o "$eur_common"             "$ct_chr" "$eur_chr" 2>> {log}
        bcftools index -f "$eur_common" 2>> {log}

        exact_ct=$(bcftools index -n "$ct_common")
        exact_eur=$(bcftools index -n "$eur_common")

        if [ "$exact_ct" -ne "$exact_eur" ]; then
            echo "ERROR: exact-match counts differ between target and reference on chr{wildcards.chrom}" >&2
            exit 1
        fi

        if bcftools query -f '%CHROM:%POS:%REF:%ALT\n' "$ct_common" 2>> {log}             | sort | uniq -d | grep -q .; then
            echo "ERROR: duplicate exact variant key in Cinque Terre common set on chr{wildcards.chrom}" >&2
            exit 1
        fi

        if bcftools query -f '%CHROM:%POS:%REF:%ALT\n' "$eur_common" 2>> {log}             | sort | uniq -d | grep -q .; then
            echo "ERROR: duplicate exact variant key in 1000G EUR common set on chr{wildcards.chrom}" >&2
            exit 1
        fi

        bcftools merge             -m none             -Ob -o {output.bcf}             "$ct_common" "$eur_common" 2>> {log}
        bcftools index -f {output.bcf} 2>> {log}

        ct_n=$(bcftools index -n "$ct_chr")
        eur_n=$(bcftools index -n "$eur_chr")
        merged_n=$(bcftools index -n {output.bcf})
        ct_samples=$(wc -l < "$tmpdir/ct.samples")
        eur_samples=$(wc -l < "$tmpdir/eur.samples")
        merged_samples=$(bcftools query -l {output.bcf} | wc -l)
        allele_mismatch=$((same_position - exact_ct))

        if [ "$merged_n" -ne "$exact_ct" ]; then
            echo "ERROR: merged variant count does not equal exact-match count on chr{wildcards.chrom}" >&2
            exit 1
        fi

        if [ "$merged_samples" -ne $((ct_samples + eur_samples)) ]; then
            echo "ERROR: merged sample count is inconsistent on chr{wildcards.chrom}" >&2
            exit 1
        fi

        {{
            printf 'chr\tCT\tEUR_QC\tsame_position\texact_match\tallele_mismatch\tCT_samples\tEUR_samples\tmerged_samples\n'
            printf 'chr{wildcards.chrom}\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'                 "$ct_n" "$eur_n" "$same_position" "$exact_ct" "$allele_mismatch"                 "$ct_samples" "$eur_samples" "$merged_samples"
        }} > {output.audit}
        """


rule concatenate_harmonized_chromosomes:
    input:
        bcfs=expand(
            "results/population_structure/harmonization/per_chromosome/chr{chrom}.harmonized.bcf",
            chrom=HARMONIZATION_CHROMS,
        ),
        indexes=expand(
            "results/population_structure/harmonization/per_chromosome/chr{chrom}.harmonized.bcf.csi",
            chrom=HARMONIZATION_CHROMS,
        )
    output:
        bcf=temp("results/population_structure/harmonization/ct_1kg_eur.harmonized.bcf"),
        index=temp("results/population_structure/harmonization/ct_1kg_eur.harmonized.bcf.csi")
    log:
        "logs/population_structure/harmonization/concat.log"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        bcftools concat             -Ob             -o {output.bcf}             {input.bcfs}             > {log} 2>&1
        bcftools index -f {output.bcf} >> {log} 2>&1
        """


rule convert_harmonized_master_to_pgen:
    input:
        bcf="results/population_structure/harmonization/ct_1kg_eur.harmonized.bcf",
        index="results/population_structure/harmonization/ct_1kg_eur.harmonized.bcf.csi"
    output:
        pgen="results/population_structure/harmonization/ct_1kg_eur.harmonized.pgen",
        pvar="results/population_structure/harmonization/ct_1kg_eur.harmonized.pvar",
        psam="results/population_structure/harmonization/ct_1kg_eur.harmonized.psam"
    log:
        "logs/population_structure/harmonization/plink2_conversion.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail

        plink2             --threads {threads}             --bcf {input.bcf}             --set-all-var-ids '@:#:$r:$a'             --make-pgen             --out results/population_structure/harmonization/ct_1kg_eur.harmonized             > {log} 2>&1
        """


rule summarize_ct_1kg_eur_harmonization:
    input:
        audits=expand(
            "results/population_structure/harmonization/audit/chr{chrom}.harmonization.tsv",
            chrom=HARMONIZATION_CHROMS,
        ),
        pvar="results/population_structure/harmonization/ct_1kg_eur.harmonized.pvar",
        psam="results/population_structure/harmonization/ct_1kg_eur.harmonized.psam",
        script="workflow/scripts/summarize_harmonization.py"
    output:
        by_chromosome="results/population_structure/harmonization/harmonization_by_chromosome.tsv",
        summary="results/population_structure/harmonization/harmonization_summary.tsv"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.script}             --audits {input.audits}             --pvar {input.pvar}             --psam {input.psam}             --by-chromosome-out {output.by_chromosome}             --summary-out {output.summary}
        """


rule write_ct_1kg_eur_harmonization_readme:
    input:
        summary="results/population_structure/harmonization/harmonization_summary.tsv"
    output:
        "results/population_structure/harmonization/README.txt"
    shell:
        r"""
        cat > {output} <<'EOF'
Cinque Terre + 1000 Genomes EUR harmonization
=============================================

Reference cohort
----------------
The EUR sample list is regenerated from the local 1000GP_Phase3.sample metadata
file using GROUP == EUR and is validated against the expected CEU, FIN, GBR,
IBS and TSI population counts.

Reference site QC
-----------------
Each Phase 3 autosomal VCF is restricted to the 503 EUR individuals and to
biallelic SNPs. AC, AN and F_MISSING are recalculated after subsetting.
Monomorphic variants, sites with >5% missing genotypes, and palindromic A/T or
C/G SNPs are removed.

The audited source callset has zero missing genotypes after these site filters,
so individual missingness is necessarily zero in the retained EUR dataset. The
workflow explicitly verifies this condition and stops if it is no longer true,
rather than silently omitting individual-level missingness QC on a changed
input.

Harmonization
-------------
The target input is the final Cinque Terre QC VCF after HWE filtering.

Only exact CHR:POS:REF:ALT matches are retained between Cinque Terre and 1000G
EUR. Same-position allele-discordant records are excluded. No strand flipping,
REF/ALT swapping, or allele-rescue procedure is used.

Large per-chromosome and concatenated BCF files are temporary Snakemake
intermediates. The persistent master harmonized dataset is:

    ct_1kg_eur.harmonized.pgen
    ct_1kg_eur.harmonized.pvar
    ct_1kg_eur.harmonized.psam

This master dataset is intentionally not MAF-filtered or LD-pruned. Those
filters belong to downstream analysis branches. PCA-specific filtering will be
applied separately.
EOF
        """


rule harmonize_1kg_eur:
    input:
        "results/reference/1kg_eur/metadata/eur_sample_inventory.tsv",
        "results/population_structure/harmonization/ct_1kg_eur.harmonized.pgen",
        "results/population_structure/harmonization/ct_1kg_eur.harmonized.pvar",
        "results/population_structure/harmonization/ct_1kg_eur.harmonized.psam",
        "results/population_structure/harmonization/harmonization_by_chromosome.tsv",
        "results/population_structure/harmonization/harmonization_summary.tsv",
        "results/population_structure/harmonization/README.txt"
