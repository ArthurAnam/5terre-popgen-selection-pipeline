# ============================================================
# Cinque Terre Population Genomics and Selection Pipeline
# ============================================================

import os
import yaml

# Version-controlled analytical configuration.
configfile: "config/config.yaml"

# Machine-specific paths are kept outside Git. When present, this file is
# merged on top of config/config.yaml so that rules never hard-code local paths.
LOCAL_CONFIG = "config/config.local.yaml"
if os.path.exists(LOCAL_CONFIG):
    with open(LOCAL_CONFIG, "r") as handle:
        local_config = yaml.safe_load(handle) or {}
    config.update(local_config)


# ============================================================
# Workflow modules
# ============================================================

include: "workflow/rules/00_preqc.smk"
include: "workflow/rules/01_sample_qc.smk"
include: "workflow/rules/99_provenance.smk"


# ============================================================
# Current reproducible target
# ============================================================
# Individual pre-QC is diagnostic at this stage: no sample is removed
# automatically from the delivered 50-sample dataset.

rule all:
    input:
        "results/preqc/preqc_done.txt",
        "results/sample_preqc/sample_preqc_summary.tsv",
        "results/provenance/software_versions.tsv"
