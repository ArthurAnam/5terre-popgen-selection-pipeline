# Pre-LASSI environment audit

The pre-LASSI environment check is intentionally split into two Snakemake
targets:

- `audit_pre_lassi_environment` generates and preserves the full diagnostic
  report. A completed audit exits successfully even when individual checks are
  recorded as `FAIL`, so Snakemake cannot delete the evidence needed to
  diagnose the failure.
- `gate_pre_lassi_environment` reads the preserved audit summary and is the
  blocking target. It exits non-zero unless
  `overall_pre_lassi_environment_status=PASS`.

Unexpected Python exceptions in the diagnostic generator still propagate as
real job failures; only *recorded analytical/environmental FAIL states* are
decoupled from report generation.

The audit covers Git/repository integrity, critical local analytical files, the
software actually visible in the active environment, SHAPEIT2 reference/map
resources, and machine resources (CPU, RAM, swap and disk). It also writes
local Conda snapshots and SHA256 fingerprints for the chromosome-level phasing
outputs.

Machine-specific reports are written under `results/provenance/` and remain
ignored by Git. After review, only portable conclusions should be copied into
`docs/analysis_records/`.

LASSI software visibility is deliberately a `WARN`, not a hard failure, at
this stage. The official LASSI package/version/scripts must be installed and
frozen only after the current machine/software/file audit has been interpreted
and closed. No production LASSI scan should depend directly on the diagnostic
target; future production scans must depend on `gate_pre_lassi_environment`.

Recommended review sequence:

1. Run `audit_pre_lassi_environment`.
2. Inspect every `FAIL` and `WARN` in
   `results/provenance/pre_lassi_environment_audit.tsv`, plus the details,
   manifest and Conda snapshots.
3. Resolve or explicitly document each issue.
4. Run `gate_pre_lassi_environment` only when the audit is expected to pass.
5. Freeze the production LASSI software, inputs and parameters only after the
   gate passes.
