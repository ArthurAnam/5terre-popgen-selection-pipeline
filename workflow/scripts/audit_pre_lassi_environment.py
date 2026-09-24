#!/usr/bin/env python3
import csv
import gzip
import hashlib
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(subprocess.check_output(["git","rev-parse","--show-toplevel"], text=True).strip())
os.chdir(ROOT)
OUT = Path("results/provenance/pre_lassi_environment_audit.tsv")
DETAIL = Path("results/provenance/pre_lassi_environment_details.txt")
MANIFEST = Path("results/provenance/pre_lassi_critical_manifest.tsv")
CONDA_LIST = Path("results/provenance/conda_list.txt")
CONDA_EXPLICIT = Path("results/provenance/conda_explicit.txt")
CONDA_HISTORY = Path("results/provenance/conda_from_history.yaml")
for p in [OUT, DETAIL, MANIFEST, CONDA_LIST, CONDA_EXPLICIT, CONDA_HISTORY]:
    p.parent.mkdir(parents=True, exist_ok=True)

rows = []
fails = 0
warns = 0
details = []

def clean(x):
    return re.sub(r"\s+", " ", str(x)).strip()

def record(check, status, value=""):
    global fails, warns
    rows.append((check, status, clean(value)))
    if status == "FAIL":
        fails += 1
    elif status == "WARN":
        warns += 1

def run(cmd, check=False):
    return subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=check)

def capture(title, cmd):
    r = run(cmd)
    details.append(f"===== {title} =====\n{r.stdout}\n")
    return r

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def gzip_ok(path):
    try:
        with gzip.open(path, "rb") as fh:
            for _ in iter(lambda: fh.read(1024 * 1024), b""):
                pass
        return True
    except Exception:
        return False

def software(name, cmd):
    exe = shutil.which(cmd[0])
    if not exe:
        record(f"software_{name}", "FAIL", f"command not found: {cmd[0]}")
        return
    r = run(cmd)
    text = " ".join(r.stdout.splitlines()[:3])
    record(f"software_{name}", "PASS" if r.returncode == 0 else "WARN", text or exe)

# Git / GitHub-local consistency
branch = run(["git","branch","--show-current"]).stdout.strip()
head = run(["git","rev-parse","HEAD"]).stdout.strip()
origin = run(["git","rev-parse","origin/main"]).stdout.strip()
record("git_branch", "PASS" if branch == "main" else "WARN", branch)
record("git_head", "INFO", head)
record("git_origin_main", "INFO", origin)
record("git_head_matches_origin", "PASS" if head == origin else "FAIL", f"HEAD={head} origin/main={origin}")
dirty = run(["git","status","--porcelain","--untracked-files=normal"]).stdout.strip()
record("git_worktree_clean", "PASS" if not dirty else "FAIL", dirty or "clean")
ab = run(["git","rev-list","--left-right","--count","origin/main...HEAD"]).stdout.strip()
record("git_ahead_behind", "INFO", ab)
fsck = capture("git fsck --full", ["git","fsck","--full"])
record("git_fsck", "PASS" if fsck.returncode == 0 else "FAIL", "repository object database OK" if fsck.returncode == 0 else fsck.stdout)
tracked = run(["git","ls-files"]).stdout.splitlines()
bad_ext = re.compile(r"\.(vcf|vcf\.gz|bcf|bed|bim|fam|haps|hap|hap\.gz|bam|cram)$", re.I)
bad = [p for p in tracked if bad_ext.search(p)]
record("git_no_large_genomic_files", "PASS" if not bad else "FAIL", "none tracked" if not bad else ";".join(bad))
for p, key in [
    ("config/config.local.yaml","git_local_config_ignored"),
    ("results/selection/phasing/shapeit2/phased/chr1.ct.maf005.phased.haps.gz","git_results_ignored"),
]:
    r = run(["git","check-ignore","-q",p])
    record(key, "PASS" if r.returncode == 0 else "FAIL", p)
capture("git status", ["git","status","--short","--branch"])
capture("git recent log", ["git","log","-n","15","--oneline","--decorate"])

# Software actually visible in the active environment
software("snakemake", ["snakemake","--version"])
software("python", ["python","--version"])
software("conda", ["conda","--version"])
software("git", ["git","--version"])
software("bcftools", ["bcftools","--version"])
software("samtools", ["samtools","--version"])
software("tabix", ["tabix","--version"])
software("bgzip", ["bgzip","--version"])
software("plink1", ["plink","--version"])
software("plink2", ["plink2","--version"])
king_exe = shutil.which("king")
if not king_exe:
    record("software_king", "FAIL", "command not found: king")
else:
    king_r = run(["king","--version"])
    king_text = " ".join(king_r.stdout.splitlines()[:3])
    king_version = re.search(r"\bKING\s+([0-9]+(?:\.[0-9]+)+)", king_text)
    if king_version:
        record("software_king", "PASS", f"{king_exe} ; KING {king_version.group(1)}")
    else:
        record(
            "software_king",
            "PASS" if king_r.returncode == 0 else "WARN",
            king_text or king_exe,
        )
software("gzip", ["gzip","--version"])

if shutil.which("smartpca"):
    eig = run(["conda","list","eigensoft"]).stdout
    ver = next((line for line in eig.splitlines() if line.startswith("eigensoft ")), "version_not_parsed")
    record("software_smartpca", "PASS", f"{shutil.which('smartpca')} ; {ver}")
else:
    record("software_smartpca", "FAIL", "smartpca not found")

try:
    import matplotlib
    record("software_matplotlib", "PASS", matplotlib.__version__)
except Exception as e:
    record("software_matplotlib", "FAIL", repr(e))

# Machine-specific phasing resources
local_cfg = Path("config/config.local.yaml")
if local_cfg.exists():
    record("local_config_exists", "PASS", str(local_cfg))
    cfg = yaml.safe_load(local_cfg.read_text()) or {}
else:
    record("local_config_exists", "FAIL", str(local_cfg))
    cfg = {}
lp = cfg.get("local_paths", {})
shapeit = Path(lp.get("shapeit2_binary","")) if lp.get("shapeit2_binary") else None
if shapeit and shapeit.is_file() and os.access(shapeit, os.X_OK):
    r = run([str(shapeit)])
    record("software_shapeit2", "PASS", " ".join(r.stdout.splitlines()[:5]))
    record("shapeit2_binary_executable", "PASS", str(shapeit))
else:
    record("software_shapeit2", "FAIL", f"configured executable missing/not executable: {shapeit}")

ref_sample = Path(lp.get("shapeit2_reference_sample","")) if lp.get("shapeit2_reference_sample") else None
record("reference_sample_file", "PASS" if ref_sample and ref_sample.is_file() else "FAIL", str(ref_sample))
group_lines = []
gp = Path("config/shapeit2_reference_groups.txt")
if gp.exists():
    group_lines = [x.strip() for x in gp.read_text().splitlines() if x.strip() and not x.lstrip().startswith("#")]
record("reference_group_scope", "PASS" if group_lines == ["EUR"] else "FAIL", ",".join(group_lines))

haps_t = lp.get("shapeit2_reference_haps_template","")
leg_t = lp.get("shapeit2_reference_legend_template","")
map_t = lp.get("shapeit2_b37_map_template","")
missing_ref = []
gzip_fail = []
for c in range(1,23):
    files = [
        Path(haps_t.format(chrom=c)) if haps_t else None,
        Path(leg_t.format(chrom=c)) if leg_t else None,
        Path(map_t.format(chrom=c)) if map_t else None,
    ]
    for f in files:
        if not f or not f.is_file() or not os.access(f, os.R_OK):
            missing_ref.append(str(f))
    for f in files[:2]:
        if f and f.is_file() and not gzip_ok(f):
            gzip_fail.append(str(f))
record("reference_files_22chr", "PASS" if not missing_ref else "FAIL", "all HAP/LEGEND/maps readable" if not missing_ref else f"{len(missing_ref)} missing/unreadable")
record("reference_gzip_integrity", "PASS" if not gzip_fail else "FAIL", "44 reference gzip files tested" if not gzip_fail else f"{len(gzip_fail)} failures")

# Critical local production files
P = Path("results/selection/phasing/shapeit2")
critical_missing = []
for c in range(1,23):
    req = [
        P/f"input/chr{c}.ct.maf005.vcf.gz",
        P/f"input/chr{c}.ct.maf005.vcf.gz.tbi",
        P/f"phased/chr{c}.ct.maf005.phased.haps.gz",
        P/f"phased/chr{c}.ct.maf005.phased.sample",
        P/f"phased/chr{c}.phase.log",
        P/f"check_after_exclude/chr{c}.check_after_exclude.ok",
    ]
    critical_missing.extend(str(f) for f in req if not f.is_file() or f.stat().st_size == 0)
record("critical_phasing_files", "PASS" if not critical_missing else "FAIL", "all 22 chromosome sets present/nonempty" if not critical_missing else f"{len(critical_missing)} missing/empty")

audit = P/"audit/phasing_audit_summary.tsv"
audit_dict = {}
if audit.exists():
    with open(audit) as fh:
        for k,v,*_ in csv.reader(fh, delimiter="\t"):
            if k != "metric":
                audit_dict[k] = v
audit_pass = (
    audit_dict.get("audit_status") == "PASS"
    and audit_dict.get("observed_total_snps") == "4914283"
    and audit_dict.get("genotype_dosage_mismatches") == "0"
    and audit_dict.get("all_logs_thread1_parsed") == "PASS"
    and audit_dict.get("lassi_structural_readiness") == "PASS"
)
record("post_phasing_audit", "PASS" if audit_pass else "FAIL", str(audit_dict) if not audit_pass else "4914283 SNP; dosage mismatch=0; all log checks PASS")

phased_gzip_fail = [str(P/f"phased/chr{c}.ct.maf005.phased.haps.gz") for c in range(1,23) if not gzip_ok(P/f"phased/chr{c}.ct.maf005.phased.haps.gz")]
record("phased_gzip_integrity", "PASS" if not phased_gzip_fail else "FAIL", "22/22 HAPS gzip OK" if not phased_gzip_fail else f"{len(phased_gzip_fail)} failures")

def bcftools_count(path):
    r = run(["bcftools","index","-n",str(path)])
    return int(r.stdout.strip()) if r.returncode == 0 and r.stdout.strip().isdigit() else None

qc_vcf = Path("results/qc/04_hwe/cinque_terre.qc_filtered.vcf.gz")
qc_n = bcftools_count(qc_vcf) if qc_vcf.exists() else None
record("final_qc_variant_count", "PASS" if qc_n == 9000246 else "FAIL", qc_n)

manifest_rows = []
maf_total = 0
for c in range(1,23):
    v = P/f"input/chr{c}.ct.maf005.vcf.gz"
    h = P/f"phased/chr{c}.ct.maf005.phased.haps.gz"
    s = P/f"phased/chr{c}.ct.maf005.phased.sample"
    l = P/f"phased/chr{c}.phase.log"
    n = bcftools_count(v)
    maf_total += n or 0
    manifest_rows.append([
        c,
        h.stat().st_size if h.exists() else 0,
        sha256(h) if h.exists() else "",
        sha256(s) if s.exists() else "",
        sha256(l) if l.exists() else "",
        n if n is not None else "",
    ])
with open(MANIFEST,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["chromosome","haps_bytes","haps_sha256","sample_sha256","phase_log_sha256","input_maf005_snps"])
    w.writerows(manifest_rows)
record("ct_maf005_input_count", "PASS" if maf_total == 5007326 else "FAIL", maf_total)
record("critical_manifest", "PASS" if all(r[2] and r[3] and r[4] for r in manifest_rows) else "FAIL", str(MANIFEST))

dry = run(["snakemake","-n","audit_selection_ct_shapeit2_phasing"])
record("phasing_audit_dag_current", "PASS" if dry.returncode == 0 else "FAIL", " ".join(dry.stdout.splitlines()[-3:]))

broken = []
for p in ROOT.rglob("*"):
    if ".git" in p.parts or ".snakemake" in p.parts:
        continue
    if p.is_symlink() and not p.exists():
        broken.append(str(p.relative_to(ROOT)))
record("broken_symlinks", "PASS" if not broken else "FAIL", "none" if not broken else ";".join(broken))

# Machine resources: informational, not arbitrary pass/fail thresholds.
mem = {}
for line in Path("/proc/meminfo").read_text().splitlines():
    if ":" in line:
        k,v=line.split(":",1)
        m=re.search(r"(\d+)",v)
        if m:
            mem[k]=int(m.group(1))*1024
record("system_cpu_logical", "INFO", os.cpu_count())
record("system_machine", "INFO", platform.platform())
record("system_memory_total_bytes", "INFO", mem.get("MemTotal",""))
record("system_memory_available_bytes", "INFO", mem.get("MemAvailable",""))
record("system_swap_total_bytes", "INFO", mem.get("SwapTotal",""))
record("system_swap_free_bytes", "INFO", mem.get("SwapFree",""))
du = shutil.disk_usage(ROOT)
record("system_repo_disk_total_bytes", "INFO", du.total)
record("system_repo_disk_free_bytes", "INFO", du.free)
record("system_load_1_5_15", "INFO", " ".join(Path("/proc/loadavg").read_text().split()[:3]))
try:
    soft,hard = __import__("resource").getrlimit(__import__("resource").RLIMIT_NOFILE)
    record("system_open_file_limit", "INFO", f"soft={soft};hard={hard}")
except Exception as e:
    record("system_open_file_limit", "WARN", repr(e))

for title,cmd in [
    ("free -h",["free","-h"]),
    ("swapon --show",["swapon","--show"]),
    ("df -hT",["df","-hT","."]),
    ("df -ih",["df","-ih","."]),
    ("lscpu",["lscpu"]),
    ("repo disk usage",["du","-sh","results","logs","benchmarks","data"]),
]:
    capture(title,cmd)

# Exact environment snapshots remain local/ignored.
snapshots = [
    (CONDA_LIST, ["conda","list"]),
    (CONDA_EXPLICIT, ["conda","list","--explicit"]),
    (CONDA_HISTORY, ["conda","env","export","--from-history"]),
]
for path,cmd in snapshots:
    r=run(cmd)
    path.write_text(r.stdout)
    record(path.stem+"_snapshot", "PASS" if r.returncode == 0 and path.stat().st_size else "WARN", str(path))

# Selected LASSI implementation: lassip v1.2.1, pinned in config/config.yaml.
lassip_cfg = lp.get("lassip_binary", "")
lassip = Path(lassip_cfg) if lassip_cfg else None
if lassip and lassip.is_file() and os.access(lassip, os.X_OK):
    r = run([str(lassip), "--help"])
    banner = " ".join(r.stdout.splitlines()[:4])
    if "lassip v1.2.1" in r.stdout:
        record("lassi_software_visibility", "PASS", f"{lassip} ; lassip v1.2.1")
    else:
        record("lassi_software_visibility", "FAIL", banner or str(lassip))
else:
    record(
        "lassi_software_visibility",
        "WARN",
        "lassip v1.2.1 selected but local binary is not yet configured; run the dedicated software preflight before production scans",
    )

record("audit_failures","INFO",fails)
record("audit_warnings","INFO",warns)
record("overall_pre_lassi_environment_status","PASS" if fails == 0 else "FAIL",f"hard_failures={fails}; warnings={warns}")

with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["check","status","value"])
    w.writerows(rows)
DETAIL.write_text("\n".join(details))

# This target is diagnostic: preserve completed reports even when checks FAIL.
# Unexpected Python exceptions still propagate as non-zero exits. The separate
# gate_pre_lassi_environment target enforces the recorded overall status.
sys.exit(0)
