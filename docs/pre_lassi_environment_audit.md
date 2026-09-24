# Pre-LASSI environment audit

The Snakemake target `audit_pre_lassi_environment` is the machine-level gate
between the completed CT phasing branch and LASSI execution.

It audits Git/repository integrity, critical local analytical files, the
software actually visible in the active environment, SHAPEIT2 reference/map
resources, and machine resources (CPU, RAM, swap and disk). It also writes
local Conda snapshots and SHA256 fingerprints for the chromosome-level phasing
outputs.

Machine-specific reports are written under `results/provenance/` and remain
ignored by Git. After review, only portable conclusions should be copied into
`docs/analysis_records/`.

LASSI software visibility is deliberately a warning rather than a hard failure
at this stage. The official LASSI script package should be installed and
frozen only after this audit confirms that the completed phasing environment
and files are internally consistent.
