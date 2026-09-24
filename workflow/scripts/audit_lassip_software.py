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

EXPECTED_REPO = "szpiech/lassip"
EXPECTED_VERSION = "1.2.1"
EXPECTED_COMMIT = "a6a9d18c2323330fbf74d5a490f9e9c4ebe41d7c"

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
        m = re.search(r'VERSION\s*=\s*"([^"]+)"', version_header.read_text())
        parsed = m.group(1) if m else ""
        record("source_version", "PASS" if parsed == EXPECTED_VERSION else "FAIL", parsed or "not parsed")
    else:
        record("source_version", "FAIL", version_header)

if binary and binary.is_file() and os.access(binary, os.X_OK):
    vr = run([str(binary), "--help"])
    m = re.search(r"lassip v([0-9.]+)", vr.stdout)
    version = m.group(1) if m else ""
    record("binary_version", "PASS" if version == EXPECTED_VERSION else "FAIL", version or vr.stdout[:300])
    record("binary_sha256", "INFO", sha256(binary))

for name, cmd in [
    ("compiler", ["g++", "--version"]),
    ("make", ["make", "--version"]),
]:
    r = run(cmd)
    record(name, "PASS" if r.returncode == 0 else "FAIL", " ".join(r.stdout.splitlines()[:2]))

record("selected_statistic", "INFO", "original LASSI T statistic via lassip --lassi")
record("selected_k", "INFO", 10)
record("selected_lassi_choice", "INFO", 4)
record("audit_failures", "INFO", fails)
record("audit_warnings", "INFO", warns)
record("overall_lassip_software_status", "PASS" if fails == 0 else "FAIL", f"hard_failures={fails}; warnings={warns}")

with OUT.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["check", "status", "value"])
    w.writerows(rows)

# Diagnostic target: always retain a completed report.
