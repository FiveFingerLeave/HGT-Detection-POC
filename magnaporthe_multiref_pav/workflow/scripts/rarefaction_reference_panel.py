#!/usr/bin/env python3
"""Section 14.1/14.2: reference-panel rarefaction.

The document specifies random sampling (n_permutations=1000) because it
is designed for a 14-genome panel, where C(14,k) is far too large to
enumerate exhaustively. Our panel was deliberately reduced to 5 host
representatives (see docs/decisions.md), so C(5,k) is at most 10 -
small enough to enumerate ALL combinations exactly instead of randomly
sampling. This is strictly more rigorous than the document's own
sampling approach, not a shortcut - documented here as a deliberate
adaptation to the smaller panel size, not a deviation in spirit.
"""
import argparse
import csv
import itertools
from statistics import mean, median


CANDIDATE_CLASSES = ["accessory_chromosome", "starship_like", "private_accessory"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--genome-ids", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-saturation", required=True)
    args = ap.parse_args()

    candidates_by_class = {c: [] for c in CANDIDATE_CLASSES}
    with open(args.manifest) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["region_type"] in CANDIDATE_CLASSES:
                carriers = set(r["reference_isolates"].split(";"))
                candidates_by_class[r["region_type"]].append(carriers)

    all_candidates = [c for cs in candidates_by_class.values() for c in cs]

    rows = []
    n = len(args.genome_ids)
    for k in range(1, n + 1):
        combos = list(itertools.combinations(args.genome_ids, k))
        for region_class, candidates in list(candidates_by_class.items()) + [("all_combined", all_candidates)]:
            counts = []
            for combo in combos:
                combo_set = set(combo)
                r_k = sum(1 for carriers in candidates if carriers & combo_set)
                counts.append(r_k)
            rows.append({
                "panel_size": k,
                "region_class": region_class,
                "mean_n_regions": f"{mean(counts):.2f}",
                "median_n_regions": f"{median(counts):.1f}",
                "ci_lower": min(counts),
                "ci_upper": max(counts),
                "n_combinations": len(combos),
            })

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # Section 14.2: saturation delta for the last panel-size step (n-1 -> n)
    sat_rows = []
    for region_class in list(candidates_by_class.keys()) + ["all_combined"]:
        r_n = next(r for r in rows if r["panel_size"] == n and r["region_class"] == region_class)
        r_n1 = next(r for r in rows if r["panel_size"] == n - 1 and r["region_class"] == region_class)
        r_n_val = float(r_n["mean_n_regions"])
        r_n1_val = float(r_n1["mean_n_regions"])
        delta = (r_n_val - r_n1_val) / r_n_val if r_n_val > 0 else 0.0
        sat_rows.append({
            "region_class": region_class,
            f"R_{n-1}_mean": r_n1_val,
            f"R_{n}_mean": r_n_val,
            "delta_fraction": f"{delta:.4f}",
            "saturated_lt_5pct": delta < 0.05,
        })

    with open(args.out_saturation, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(sat_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(sat_rows)

    print(f"Wrote {len(rows)} rarefaction rows to {args.out}")
    print(f"Wrote saturation summary to {args.out_saturation}")
    for r in sat_rows:
        print(f"  {r['region_class']}: delta={r['delta_fraction']} saturated={r['saturated_lt_5pct']}")


if __name__ == "__main__":
    main()
