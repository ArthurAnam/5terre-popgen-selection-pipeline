# ============================================================
# Runtime software provenance
# ============================================================
# Records the versions actually visible in the environment used to run the
# workflow. This complements envs/pipeline.yaml, which describes dependencies.

rule record_software_versions:
    output:
        "results/provenance/software_versions.tsv"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/provenance

        {
            printf "software\tversion\n"
            printf "snakemake\t%s\n" "$(snakemake --version | head -n 1)"
            printf "python\t%s\n" "$(python --version 2>&1 | cut -d ' ' -f 2)"
            printf "bcftools\t%s\n" "$(bcftools --version | head -n 1 | cut -d ' ' -f 2)"
            printf "htslib\t%s\n" "$(bcftools --version | sed -n 's/^Using htslib //p' | head -n 1)"
            printf "tabix\t%s\n" "$(tabix --version 2>&1 | head -n 1 | grep -oE '[0-9]+([.][0-9]+)+' | head -n 1)"
            printf "git\t%s\n" "$(git --version | cut -d ' ' -f 3)"
            printf "conda\t%s\n" "$(conda --version | cut -d ' ' -f 2)"
        } > {output}
        """
