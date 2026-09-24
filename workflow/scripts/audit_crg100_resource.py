#!/usr/bin/env python3
import csv
import hashlib
import os
import re
from pathlib import Path

import yaml

OUT = Path("results/selection/lassi/preflight/crg100_resource_audit.tsv")
OUT.parent.mkdir(parents=True, exist_ok=True)

EXPECTED_FILENAME = "wgEncodeCrgMapabilityAlign100mer.bigWig"
EXPECTED_MD5 = "a1b1a8c99431fedf6a3b4baef028cca4"
EXPECTED_BUILD = "GRCh37/hg19"
EXPECTED_UCSC_ACCESSION = "wgEncodeEH000317"
MIN_BYTES = 90 * 1024 * 1024

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

def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

main_cfg_path = Path("config/config.yaml")
local_cfg_path = Path("config/config.local.yaml")

if not main_cfg_path.is_file():
    record("main_config", "FAIL", main_cfg_path)
    mcfg = {}
else:
    record("main_config", "PASS", main_cfg_path)
    mcfg = yaml.safe_load(main_cfg_path.read_text()) or {}

mapp = (((mcfg.get("selection", {}) or {}).get("lassi", {}) or {}).get("mappability", {}) or {})
record("configured_genome_build", "PASS" if mapp.get("genome_build") == EXPECTED_BUILD else "FAIL", mapp.get("genome_build", "missing"))
record("configured_expected_filename", "PASS" if mapp.get("expected_filename") == EXPECTED_FILENAME else "FAIL", mapp.get("expected_filename", "missing"))
record("configured_expected_md5", "PASS" if mapp.get("expected_md5") == EXPECTED_MD5 else "FAIL", mapp.get("expected_md5", "missing"))
record("configured_ucsc_accession", "PASS" if mapp.get("ucsc_accession") == EXPECTED_UCSC_ACCESSION else "FAIL", mapp.get("ucsc_accession", "missing"))
record("configured_threshold", "PASS" if float(mapp.get("threshold", -1)) == 0.9 else "FAIL", mapp.get("threshold", "missing"))

if not local_cfg_path.is_file():
    record("local_config", "FAIL", local_cfg_path)
    local_paths = {}
else:
    record("local_config", "PASS", local_cfg_path)
    lcfg = yaml.safe_load(local_cfg_path.read_text()) or {}
    local_paths = lcfg.get("local_paths", {}) or {}

raw_path = local_paths.get("crg100_bigwig")
path = Path(raw_path).expanduser() if raw_path else None
record("crg100_path_configured", "PASS" if path else "FAIL", raw_path or "missing")
record("crg100_exists", "PASS" if path and path.is_file() else "FAIL", path or "missing")

if path and path.is_file():
    record("crg100_filename", "PASS" if path.name == EXPECTED_FILENAME else "WARN", path.name)
    size = path.stat().st_size
    record("crg100_size_bytes", "PASS" if size >= MIN_BYTES else "FAIL", size)
    observed_md5 = md5sum(path)
    record("crg100_md5", "PASS" if observed_md5 == EXPECTED_MD5 else "FAIL", observed_md5)
    with open(path, "rb") as fh:
        magic = fh.read(4).hex()
    # BigWig magic bytes may appear little-endian as 26fc8f88.
    record("crg100_bigwig_magic", "PASS" if magic in {"26fc8f88", "888ffc26"} else "FAIL", magic)

record("audit_failures", "INFO", fails)
record("audit_warnings", "INFO", warns)
record("overall_crg100_resource_status", "PASS" if fails == 0 else "FAIL", f"hard_failures={fails}; warnings={warns}")

with OUT.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["check", "status", "value"])
    w.writerows(rows)
