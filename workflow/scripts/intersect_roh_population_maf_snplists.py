#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def read_ids(path):
    with open(path, "r", encoding="utf-8") as h:
        return {line.strip() for line in h if line.strip()}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--joint-snplist", required=True)
    p.add_argument("--population-list", action="append", required=True,
                   help="POP=path/to/snplist; repeat once per population")
    p.add_argument("--intersection-out", required=True)
    p.add_argument("--counts-out", required=True)
    p.add_argument("--summary-out", required=True)
    args = p.parse_args()

    joint_order = []
    with open(args.joint_snplist, "r", encoding="utf-8") as h:
        joint_order = [line.strip() for line in h if line.strip()]
    joint_set = set(joint_order)

    pop_sets = {}
    for spec in args.population_list:
        pop, path = spec.split("=", 1)
        ids = read_ids(path)
        if not ids.issubset(joint_set):
            raise SystemExit(f"ERROR: {pop} MAF list contains IDs outside joint MAF panel")
        pop_sets[pop] = ids

    intersection = set.intersection(*(pop_sets[p] for p in sorted(pop_sets)))
    ordered_intersection = [vid for vid in joint_order if vid in intersection]

    Path(args.intersection_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.intersection_out, "w", encoding="utf-8") as out:
        for vid in ordered_intersection:
            out.write(vid + "\n")

    with open(args.counts_out, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["population", "n_variants_maf_ge_0.05", "fraction_of_joint_maf_panel"])
        for pop in sorted(pop_sets):
            w.writerow([pop, len(pop_sets[pop]), len(pop_sets[pop]) / len(joint_order)])

    with open(args.summary_out, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["joint_maf0.05_variants", len(joint_order)])
        w.writerow(["populations", ",".join(sorted(pop_sets))])
        w.writerow(["all_populations_maf0.05_intersection_variants", len(ordered_intersection)])
        w.writerow(["intersection_fraction_of_joint_maf_panel", len(ordered_intersection) / len(joint_order)])

if __name__ == "__main__":
    main()
