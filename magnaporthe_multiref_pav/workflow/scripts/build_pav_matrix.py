#!/usr/bin/env python3
"""Combine per-sample region-level PAV calls (pav_call.py output) into
one panel_id x sample_id matrix."""
import argparse
import csv
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region-calls", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    matrix = defaultdict(dict)
    region_type = {}
    samples = []

    for path in args.region_calls:
        with open(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                sid = r["sample_id"]
                if sid not in samples:
                    samples.append(sid)
                matrix[r["panel_id"]][sid] = r["call"]
                region_type[r["panel_id"]] = r["region_type"]

    with open(args.out, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["panel_id", "region_type"] + samples)
        for panel_id in sorted(matrix.keys()):
            row = [panel_id, region_type[panel_id]]
            row += [matrix[panel_id].get(s, "NA") for s in samples]
            writer.writerow(row)

    print(f"Wrote PAV matrix: {len(matrix)} panel regions x {len(samples)} samples -> {args.out}")


if __name__ == "__main__":
    main()
