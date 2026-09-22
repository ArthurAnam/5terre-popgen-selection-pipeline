#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def read_annotations(path):
    mapping = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "population"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit("ERROR: annotations must contain sample_id and population")
        for row in reader:
            sid = row["sample_id"]
            if sid in mapping:
                raise SystemExit(f"ERROR: duplicate annotation for {sid}")
            mapping[sid] = row["population"]
    return mapping


def read_psam(path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        header = None
        iid_idx = None
        fid_idx = None
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.split()
            if header is None:
                header = [x.lstrip("#") for x in fields]
                if "IID" not in header:
                    raise SystemExit("ERROR: PSAM header lacks IID")
                iid_idx = header.index("IID")
                fid_idx = header.index("FID") if "FID" in header else None
                continue
            iid = fields[iid_idx]
            fid = fields[fid_idx] if fid_idx is not None else "0"
            rows.append((fid, iid))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psam", required=True)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    annotations = read_annotations(args.annotations)
    rows = read_psam(args.psam)

    seen = set()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        for fid, iid in rows:
            if iid in seen:
                raise SystemExit(f"ERROR: duplicate IID in PSAM: {iid}")
            seen.add(iid)
            if iid not in annotations:
                raise SystemExit(f"ERROR: missing population annotation for {iid}")
            handle.write(f"{fid}\t{iid}\t{annotations[iid]}\n")

    missing = sorted(set(annotations) - seen)
    if missing:
        raise SystemExit(f"ERROR: {len(missing)} annotated samples absent from PSAM")


if __name__ == "__main__":
    main()
