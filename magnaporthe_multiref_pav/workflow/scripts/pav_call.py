#!/usr/bin/env python3
"""Section 11: window-based PAV calling for one pilot isolate against
the panel, combining two evaluation tiers (Section 10.2):

- "unique" tier: MAPQ-filtered, primary-only mosdepth (secondary/
  supplementary excluded) - the quantitative PAV evidence.
- "all" tier: unfiltered, secondary/supplementary INCLUDED - used only
  to flag windows whose apparent presence collapses once secondaries
  are excluded (repeat_ambiguous / multi-mapping artefact warning,
  exactly the check the document requires before trusting a
  high-homology region's primary-alignment call alone).

Panel FASTA headers ARE the mosdepth "chrom" field (no whitespace in
our header convention), so panel_id/region_type are parsed directly
from each window's chrom - no separate manifest join needed.
"""
import argparse
import csv
import gzip
from collections import defaultdict


def parse_chrom(chrom):
    # >PANEL000001|type=strict_core|cluster=CORE_001|rep=7015|source=...
    parts = dict(p.split("=", 1) for p in chrom.split("|")[1:] if "=" in p)
    panel_id = chrom.split("|", 1)[0]
    return panel_id, parts.get("type", "unclassified")


def load_thresholds(path, depth_cutoff_index=0):
    """mosdepth --thresholds output: chrom start end name thresh1_bases[,thresh2_bases,...]
    Returns dict[(chrom,start,end)] -> bases_at_or_above_cutoff."""
    out = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        header = f.readline()  # mosdepth writes a header line
        for line in f:
            cols = line.rstrip("\n").split("\t")
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            bases = int(cols[4 + depth_cutoff_index])
            out[(chrom, start, end)] = bases
    return out


def load_mean_depth(path):
    out = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            mean_depth = float(cols[-1])
            out[(chrom, start, end)] = mean_depth
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--unique-thresholds", required=True)
    ap.add_argument("--unique-regions", required=True)
    ap.add_argument("--all-thresholds", required=True)
    ap.add_argument("--present-breadth-core", type=float, required=True)
    ap.add_argument("--present-breadth-accessory", type=float, required=True)
    ap.add_argument("--absent-breadth", type=float, required=True)
    ap.add_argument("--ambiguous-gap", type=float, default=0.3)
    ap.add_argument("--out-windows", required=True)
    ap.add_argument("--out-regions", required=True)
    args = ap.parse_args()

    core_types = {"strict_core", "soft_core"}

    unique_thresh = load_thresholds(args.unique_thresholds)
    all_thresh = load_thresholds(args.all_thresholds)
    unique_depth = load_mean_depth(args.unique_regions)

    window_rows = []
    region_windows = defaultdict(list)

    for key, bases_unique in unique_thresh.items():
        chrom, start, end = key
        panel_id, region_type = parse_chrom(chrom)
        length = end - start
        breadth_unique = bases_unique / length if length else 0.0
        bases_all = all_thresh.get(key, bases_unique)
        breadth_all = bases_all / length if length else 0.0
        mean_depth = unique_depth.get(key, 0.0)

        present_cut = args.present_breadth_core if region_type in core_types else args.present_breadth_accessory
        multimapping_gap = breadth_all - breadth_unique

        if multimapping_gap > args.ambiguous_gap:
            call = "ambiguous_multimapping"
        elif breadth_unique >= present_cut:
            call = "present"
        elif breadth_unique < args.absent_breadth:
            call = "absent"
        else:
            call = "uncertain"

        window_rows.append({
            "sample_id": args.sample_id, "panel_id": panel_id, "region_type": region_type,
            "chrom": chrom, "start": start, "end": end,
            "mean_depth_unique": f"{mean_depth:.2f}",
            "breadth_unique": f"{breadth_unique:.4f}", "breadth_all": f"{breadth_all:.4f}",
            "call": call,
        })
        region_windows[(panel_id, region_type)].append(call)

    with open(args.out_windows, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(window_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(window_rows)

    region_rows = []
    for (panel_id, region_type), calls in region_windows.items():
        n = len(calls)
        counts = {c: calls.count(c) for c in set(calls)}
        majority_call = max(counts, key=counts.get)
        region_rows.append({
            "sample_id": args.sample_id, "panel_id": panel_id, "region_type": region_type,
            "n_windows": n, "call": majority_call,
            "n_present": counts.get("present", 0), "n_absent": counts.get("absent", 0),
            "n_uncertain": counts.get("uncertain", 0),
            "n_ambiguous_multimapping": counts.get("ambiguous_multimapping", 0),
        })

    with open(args.out_regions, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(region_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(region_rows)

    print(f"{args.sample_id}: {len(window_rows)} windows, {len(region_rows)} panel regions")


if __name__ == "__main__":
    main()
