rule all:
    input:
        "results/dummy.txt"


rule dummy:
    output:
        "results/dummy.txt"
    shell:
        "echo 'Pipeline works' > {output}"
