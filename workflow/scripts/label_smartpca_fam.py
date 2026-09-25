#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def read_annotations(path: Path):
    mapping = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sid = row["sample_id"]
            if sid in mapping:
                raise SystemExit(f"ERROR: duplicate annotation for sample {sid}")
            mapping[sid] = row["population"]
    return mapping


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fam", required=True)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--output-fam", required=True)
    args = parser.parse_args()

    annotations = read_annotations(Path(args.annotations))
    seen = set()
    output = Path(args.output_fam)
    output.parent.mkdir(parents=True, exist_ok=True)

    with Path(args.fam).open("r", encoding="utf-8") as src, output.open(
        "w", encoding="utf-8"
    ) as dst:
        for raw in src:
            if not raw.strip():
                continue
            fields = raw.split()
            if len(fields) < 6:
                raise SystemExit("ERROR: malformed PLINK FAM row")
            iid = fields[1]
            if iid not in annotations:
                raise SystemExit(f"ERROR: no population annotation for FAM sample {iid}")
            if iid in seen:
                raise SystemExit(f"ERROR: duplicate FAM sample ID {iid}")
            seen.add(iid)
            fields[5] = annotations[iid]
            dst.write("\t".join(fields[:6]) + "\n")

    missing = sorted(set(annotations) - seen)
    if missing:
        raise SystemExit(
            f"ERROR: {len(missing)} annotated samples are absent from the PLINK FAM"
        )


if __name__ == "__main__":
    main()
