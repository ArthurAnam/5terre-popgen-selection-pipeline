#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path

def parse_expected(text):
    out = {}
    for item in text.split(","):
        pop, n = item.split(":", 1)
        out[pop] = int(n)
    return out

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fam", required=True)
    p.add_argument("--annotations", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--expected-counts", required=True)
    args = p.parse_args()

    expected = parse_expected(args.expected_counts)
    annotation = {}
    with open(args.annotations, "r", encoding="utf-8") as h:
        r = csv.DictReader(h, delimiter="\t")
        if not {"sample_id", "population"}.issubset(r.fieldnames or []):
            raise SystemExit("ERROR: annotations require sample_id and population")
        for row in r:
            annotation[row["sample_id"]] = row["population"]

    grouped = defaultdict(list)
    with open(args.fam, "r", encoding="utf-8") as h:
        for raw in h:
            if not raw.strip():
                continue
            fields = raw.split()
            fid, iid = fields[0], fields[1]
            if iid not in annotation:
                raise SystemExit(f"ERROR: missing population annotation for {iid}")
            grouped[annotation[iid]].append((fid, iid))

    observed_pops = set(grouped)
    if observed_pops != set(expected):
        raise SystemExit(f"ERROR: population set mismatch: observed={sorted(observed_pops)} expected={sorted(expected)}")

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    for pop in sorted(expected):
        rows = grouped[pop]
        if len(rows) != expected[pop]:
            raise SystemExit(f"ERROR: {pop}: expected {expected[pop]} samples, found {len(rows)}")
        with open(outdir / f"{pop}.keep", "w", encoding="utf-8") as out:
            for fid, iid in rows:
                out.write(f"{fid}\t{iid}\n")

if __name__ == "__main__":
    main()
