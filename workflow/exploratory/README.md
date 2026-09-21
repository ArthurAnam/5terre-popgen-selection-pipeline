# Exploratory QC scripts

These scripts reproduce diagnostic analyses used during development to inspect
the input data and justify QC decisions.

They are intentionally separated from `workflow/scripts/`, which contains code
required by the default production analysis. Their Snakemake rules remain
available as explicit targets when the diagnostic history needs to be
reproduced, audited, or revisited for reviewer questions.

They are not part of the default `rule all`.
