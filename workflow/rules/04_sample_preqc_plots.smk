rule sample_preqc_plots:
    input:
        metrics="results/sample_preqc/sample_metrics_autosomal.tsv",
        genotype_fields="results/sample_preqc/sample_genotype_field_diagnostics.tsv"
    output:
        depth_missing="results/sample_preqc/plots/depth_vs_missingness.png",
        het_ab="results/sample_preqc/plots/heterozygosity_vs_allelic_imbalance.png",
        readme="results/sample_preqc/plots/README.txt"
    conda:
        "../../envs/pipeline.yaml"
    shell:
        r"""
        set -euo pipefail
        python workflow/exploratory/plot_sample_preqc.py
        """
