# ============================================================
# Pairwise population differentiation
# ============================================================
# Hudson FST is computed on the same high-LD-masked, MAF-filtered,
# LD-pruned joint panel used for the PCA sensitivity analysis.
# The former average-linkage clustering branch is intentionally not retained.

FST_DIR = "results/population_structure/fst"
FST_WITHIN = FST_DIR + "/population_clusters.txt"
FST_PREFIX = FST_DIR + "/pairwise_hudson"


rule prepare_pairwise_fst_clusters:
    input:
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        annotations=PCA_SAMPLE_ANNOTATIONS
    output:
        clusters=FST_WITHIN
    params:
        expected_samples=lambda wc: config["population_structure"]["fst"]["expected_samples"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {FST_DIR}

        awk '
        BEGIN {{FS=OFS="\t"}}
        FNR==NR {{
            if (FNR>1) pop[$1]=$2
            next
        }}
        FNR==1 {{
            for (i=1; i<=NF; i++) {{
                h=$i
                sub(/^#/, "", h)
                if (h=="FID") fid_col=i
                if (h=="IID") iid_col=i
            }}
            if (!iid_col) {{
                print "ERROR: PSAM header lacks IID" > "/dev/stderr"
                exit 2
            }}
            next
        }}
        {{
            iid=$iid_col
            fid=(fid_col ? $fid_col : "0")
            if (!(iid in pop)) {{
                print "ERROR: missing population annotation for " iid > "/dev/stderr"
                exit 3
            }}
            print fid, iid, pop[iid]
            seen[iid]=1
        }}
        END {{
            for (iid in pop) {{
                if (!(iid in seen)) {{
                    print "ERROR: annotated sample absent from PSAM: " iid > "/dev/stderr"
                    bad=1
                }}
            }}
            if (bad) exit 4
        }}
        ' {input.annotations} {input.psam} > {output.clusters}

        n=$(wc -l < {output.clusters})
        [ "$n" -eq "{params.expected_samples}" ] || {{
            echo "ERROR: expected {params.expected_samples} FST samples, found $n" >&2
            exit 1
        }}
        """


rule pairwise_hudson_fst:
    input:
        pgen=PCA_HIGHLD_PRUNED_PREFIX + ".pgen",
        pvar=PCA_HIGHLD_PRUNED_PREFIX + ".pvar",
        psam=PCA_HIGHLD_PRUNED_PREFIX + ".psam",
        clusters=FST_WITHIN
    output:
        summary=FST_PREFIX + ".fst.summary"
    params:
        blocksize=lambda wc: config["population_structure"]["fst"]["blocksize_variants"],
        expected_variants=lambda wc: config["population_structure"]["fst"]["expected_variants"]
    log:
        "logs/population_structure/fst/pairwise_hudson.log"
    threads: 1
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {FST_DIR} logs/population_structure/fst

        n_variants=$(awk 'BEGIN{{n=0}} !/^#/ {{n++}} END{{print n}}' {input.pvar})
        [ "$n_variants" -eq "{params.expected_variants}" ] || {{
            echo "ERROR: expected {params.expected_variants} FST markers, found $n_variants" >&2
            exit 1
        }}

        plink2 \
            --threads {threads} \
            --pfile {PCA_HIGHLD_PRUNED_PREFIX} \
            --within {input.clusters} PCA_POP \
            --fst PCA_POP method=hudson blocksize={params.blocksize} cols=+nobs \
            --out {FST_PREFIX} \
            > {log} 2>&1

        test -s {output.summary}
        n_pairs=$(awk 'NR>1 {{n++}} END{{print n+0}}' {output.summary})
        [ "$n_pairs" -eq 15 ] || {{
            echo "ERROR: expected 15 pairwise FST estimates, found $n_pairs" >&2
            exit 1
        }}
        """


rule run_pairwise_fst:
    input:
        FST_PREFIX + ".fst.summary"
