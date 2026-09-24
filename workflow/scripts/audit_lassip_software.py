#!/usr/bin/env python3
import csv
import hashlib
import os
import re
import subprocess
from pathlib import Path

import yaml

ROOT = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
os.chdir(ROOT)
OUT = Path("results/selection/lassi/preflight/lassip_software_audit.tsv")
OUT.parent.mkdir(parents=True, exist_ok=True)

# Frozen production expectations. Intentional method changes should update both
# config/config.yaml and this audit contract.
EXPECTED_REPO = "szpiech/lassip"
EXPECTED_VERSION = "1.2.1"
EXPECTED_COMMIT = "a6a9d18c2323330fbf74d5a490f9e9c4ebe41d7c"
EXPECTED_PRIMARY_METHOD = "saltiLASSI"
EXPECTED_PRIMARY_STATISTIC = "Lambda"
EXPECTED_PRIMARY_FLAG = "--salti"
EXPECTED_SECONDARY_FLAG = "--lassi"
EXPECTED_K = 10
EXPECTED_LASSI_CHOICE = 4
EXPECTED_DIST_TYPE = "bp"
EXPECTED_MAX_EXTEND_BP = 100000

rows = []
fails = 0
warns = 0

def record(check, status, value=""):
    global fails, warns
    rows.append((check, status, re.sub(r"\s+", " ", str(value)).strip()))
    if status == "FAIL":
        fails += 1
    elif status == "WARN":
        warns += 1

def run(cmd):
    return subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

# Tracked scientific configuration.
main_cfg_path = Path("config/config.yaml")
if not main_cfg_path.is_file():
    record("main_config", "FAIL", main_cfg_path)
    lassi_cfg = {}
else:
    main_cfg = yaml.safe_load(main_cfg_path.read_text()) or {}
    lassi_cfg = (main_cfg.get("selection", {}) or {}).get("lassi", {}) or {}
    record("main_config", "PASS", main_cfg_path)

impl = lassi_cfg.get("implementation", {}) or {}
primary = lassi_cfg.get("primary_model", {}) or {}
secondary = lassi_cfg.get("secondary_lassi_T", {}) or {}

record(
    "configured_primary_method",
    "PASS" if lassi_cfg.get("method") == EXPECTED_PRIMARY_METHOD else "FAIL",
    lassi_cfg.get("method", "missing"),
)
record(
    "configured_primary_statistic",
    "PASS" if lassi_cfg.get("statistic") == EXPECTED_PRIMARY_STATISTIC else "FAIL",
    lassi_cfg.get("statistic", "missing"),
)
record(
    "configured_primary_mode_flag",
    "PASS" if impl.get("primary_mode_flag") == EXPECTED_PRIMARY_FLAG else "FAIL",
    impl.get("primary_mode_flag", "missing"),
)
record(
    "configured_secondary_mode_flag",
    "PASS" if secondary.get("enabled") is True and impl.get("secondary_mode_flag") == EXPECTED_SECONDARY_FLAG else "FAIL",
    f"enabled={secondary.get('enabled')}; flag={impl.get('secondary_mode_flag', 'missing')}",
)
record(
    "configured_k",
    "PASS" if primary.get("k_truncation") == EXPECTED_K else "FAIL",
    primary.get("k_truncation", "missing"),
)
record(
    "configured_lassi_choice",
    "PASS" if primary.get("lassi_choice") == EXPECTED_LASSI_CHOICE else "FAIL",
    primary.get("lassi_choice", "missing"),
)
record(
    "configured_distance_type",
    "PASS" if primary.get("distance_type") == EXPECTED_DIST_TYPE else "FAIL",
    primary.get("distance_type", "missing"),
)
record(
    "configured_max_extend_bp",
    "PASS" if primary.get("max_extend_bp") == EXPECTED_MAX_EXTEND_BP else "FAIL",
    primary.get("max_extend_bp", "missing"),
)
record(
    "configured_software_version",
    "PASS" if impl.get("version") == EXPECTED_VERSION else "FAIL",
    impl.get("version", "missing"),
)
record(
    "configured_software_commit",
    "PASS" if impl.get("commit") == EXPECTED_COMMIT else "FAIL",
    impl.get("commit", "missing"),
)

# Machine-local paths.
cfg_path = Path("config/config.local.yaml")
if not cfg_path.is_file():
    record("local_config", "FAIL", cfg_path)
    lp = {}
else:
    record("local_config", "PASS", cfg_path)
    cfg = yaml.safe_load(cfg_path.read_text()) or {}
    lp = cfg.get("local_paths", {})

src = Path(lp.get("lassip_source_dir", "")) if lp.get("lassip_source_dir") else None
binary = Path(lp.get("lassip_binary", "")) if lp.get("lassip_binary") else None

record("source_dir", "PASS" if src and src.is_dir() else "FAIL", src)
record("binary_exists", "PASS" if binary and binary.is_file() else "FAIL", binary)
record("binary_executable", "PASS" if binary and binary.is_file() and os.access(binary, os.X_OK) else "FAIL", binary)

if src and src.is_dir():
    origin = run(["git", "-C", str(src), "remote", "get-url", "origin"])
    origin_text = origin.stdout.strip()
    origin_ok = origin.returncode == 0 and (
        "github.com/szpiech/lassip" in origin_text or origin_text.endswith("szpiech/lassip.git")
    )
    record("upstream_repository", "PASS" if origin_ok else "FAIL", origin_text)

    head = run(["git", "-C", str(src), "rev-parse", "HEAD"])
    head_text = head.stdout.strip()
    record("upstream_commit", "PASS" if head.returncode == 0 and head_text == EXPECTED_COMMIT else "FAIL", head_text)

    dirty = run(["git", "-C", str(src), "status", "--porcelain", "--untracked-files=no"])
    record("tracked_source_clean", "PASS" if dirty.returncode == 0 and not dirty.stdout.strip() else "FAIL", dirty.stdout.strip() or "clean")

    version_header = src / "src" / "lassip-cli.h"
    if version_header.is_file():
        header_text = version_header.read_text()
        m = re.search(r'VERSION\s*=\s*"([^"]+)"', header_text)
        parsed = m.group(1) if m else ""
        record("source_version", "PASS" if parsed == EXPECTED_VERSION else "FAIL", parsed or "not parsed")
        default_bp = re.search(r"DEFAULT_MAX_EXTEND_BP\s*=\s*([0-9.]+)", header_text)
        parsed_bp = float(default_bp.group(1)) if default_bp else None
        record(
            "source_default_max_extend_bp",
            "PASS" if parsed_bp == float(EXPECTED_MAX_EXTEND_BP) else "FAIL",
            parsed_bp if parsed_bp is not None else "not parsed",
        )
    else:
        record("source_version", "FAIL", version_header)
        record("source_default_max_extend_bp", "FAIL", version_header)

if binary and binary.is_file() and os.access(binary, os.X_OK):
    vr = run([str(binary), "--help"])
    help_text = vr.stdout
    m = re.search(r"lassip v([0-9.]+)", help_text)
    version = m.group(1) if m else ""
    record("binary_version", "PASS" if version == EXPECTED_VERSION else "FAIL", version or help_text[:300])
    record("binary_sha256", "INFO", sha256(binary))
    record("binary_supports_salti", "PASS" if "--salti" in help_text else "FAIL", "--salti")
    record("binary_supports_lassi", "PASS" if "--lassi" in help_text else "FAIL", "--lassi")
    record("binary_supports_dist_type", "PASS" if "--dist-type" in help_text else "FAIL", "--dist-type")
    record("binary_supports_max_extend_bp", "PASS" if "--max-extend-bp" in help_text else "FAIL", "--max-extend-bp")
else:
    record("binary_version", "FAIL", "binary unavailable")
    record("binary_supports_salti", "FAIL", "--salti")
    record("binary_supports_lassi", "FAIL", "--lassi")
    record("binary_supports_dist_type", "FAIL", "--dist-type")
    record("binary_supports_max_extend_bp", "FAIL", "--max-extend-bp")

for name, cmd in [
    ("compiler", ["g++", "--version"]),
    ("make", ["make", "--version"]),
]:
    r = run(cmd)
    record(name, "PASS" if r.returncode == 0 else "FAIL", " ".join(r.stdout.splitlines()[:2]))

record("selected_statistic", "PASS", "saltiLASSI Lambda via lassip --salti")
record("secondary_statistic", "PASS", "original LASSI T via lassip --lassi")
record("selected_k", "PASS", EXPECTED_K)
record("selected_lassi_choice", "PASS", EXPECTED_LASSI_CHOICE)
record("selected_distance_type", "PASS", EXPECTED_DIST_TYPE)
record("selected_max_extend_bp", "PASS", EXPECTED_MAX_EXTEND_BP)
record("audit_failures", "INFO", fails)
record("audit_warnings", "INFO", warns)
record("overall_lassip_software_status", "PASS" if fails == 0 else "FAIL", f"hard_failures={fails}; warnings={warns}")

with OUT.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["check", "status", "value"])
    w.writerows(rows)

# Diagnostic target: always retain a completed report.
