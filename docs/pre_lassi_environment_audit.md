# Pre-LASSI environment and software gates

The pre-LASSI checks are split into diagnostic audits and blocking gates so that failed diagnostics remain available for inspection rather than being deleted as failed Snakemake outputs.

## Environment audit

`audit_pre_lassi_environment` records repository state, active software, critical local phasing/reference files, gzip integrity, machine resources and Conda snapshots.

`gate_pre_lassi_environment` is the blocking target and passes only when the recorded environment conditions are acceptable.

Machine-specific reports live under `results/provenance/` and are not committed. Portable conclusions are copied into `docs/analysis_records/`.

The completed environment audit/gate is part of provenance for the current phased CT dataset.

## lassip software gate

The earlier environment audit intentionally treated legacy `LASSI_iterator.py` visibility as a warning because the production implementation had not yet been frozen. That warning is historical and must not be interpreted as the current software state.

Production software is now frozen separately as:

* `lassip` v1.2.1
* commit `a6a9d18c2323330fbf74d5a490f9e9c4ebe41d7c`
* saltiLASSI Lambda primary via `--salti`
* original LASSI T secondary via `--lassi`

The dedicated lassip software gate has passed and its binary SHA256 is recorded in `docs/analysis_records/selection_saltilassi_method_freeze.tsv`.

## CRG100 resource gate

The production hg19 CRG100 BigWig is also handled by a dedicated audit/gate. The resource gate verifies the expected file, checksum and BigWig signature before downstream use.

No production scan should bypass these gates.
