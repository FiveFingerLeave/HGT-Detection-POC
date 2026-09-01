#!/usr/bin/env python3
"""Presence/absence (PAV) calling from mosdepth coverage output, per the
POC guideline (Dokumentation/POC Workflow ... .md, Section 1).

A region is called:
- "present"   if breadth >= --present-breadth AND mean depth >= --present-depth
- "absent"    if breadth <= --absent-breadth
- "uncertain" otherwise (the gray zone between the two breadth thresholds,
  or high breadth with insufficient depth)

Breadth is bases-covered-at-or-above-`--present-depth` divided by region
length, taken from mosdepth's `--thresholds <present_depth>` output (the
`.thresholds.bed.gz` file), not from `.regions.bed.gz` alone (which only
gives mean depth, not the fraction of the region actually covered - the
two are not interchangeable: a region can have high mean depth from a
small highly-covered sub-region while mostly being uncovered).
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def read_regions(path: Path) -> dict[str, dict[str, object]]:
    """Return {region_id: {"length_bp": int, "region_class": str}} from a
    BED file (chrom, start, end, region_id[, region_class]). region_class
    (5th column) distinguishes "core" marker loci from "candidate"
    HGT/mChr/Starship loci (Section 3 of the POC guideline); it is "" if
    the BED has no 5th column.
    """
    regions: dict[str, dict[str, object]] = {}
    with path.open("r") as infile:
        for line in infile:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end, region_id = fields[0], int(fields[1]), int(fields[2]), fields[3]
            region_class = fields[4] if len(fields) > 4 else ""
            regions[region_id] = {"length_bp": end - start, "region_class": region_class}
    return regions


def _open_maybe_gzip(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r")


def read_mean_depth(path: Path) -> dict[str, float]:
    """Return {region_id: mean_depth} from mosdepth's <prefix>.regions.bed.gz."""
    depths: dict[str, float] = {}
    with _open_maybe_gzip(path) as infile:
        for line in infile:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            region_id, mean_depth = fields[3], float(fields[4])
            depths[region_id] = mean_depth
    return depths


def read_threshold_bases(path: Path) -> dict[str, int]:
    """Return {region_id: bases_at_or_above_threshold} from mosdepth's
    <prefix>.thresholds.bed.gz (single-threshold case: one count column)."""
    bases: dict[str, int] = {}
    with _open_maybe_gzip(path) as infile:
        for line in infile:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            region_id, count = fields[3], int(fields[4])
            bases[region_id] = count
    return bases


def call_pav(
    sample_id: str,
    regions: dict[str, dict[str, object]],
    mean_depths: dict[str, float],
    threshold_bases: dict[str, int],
    present_breadth: float,
    present_depth: float,
    absent_breadth: float,
) -> list[dict[str, object]]:
    rows = []
    for region_id, region_info in regions.items():
        length_bp = region_info["length_bp"]
        mean_depth = mean_depths.get(region_id, 0.0)
        bases_covered = threshold_bases.get(region_id, 0)
        breadth = bases_covered / length_bp if length_bp > 0 else 0.0

        if breadth >= present_breadth and mean_depth >= present_depth:
            call = "present"
        elif breadth <= absent_breadth:
            call = "absent"
        else:
            call = "uncertain"

        rows.append(
            {
                "sample_id": sample_id,
                "region_id": region_id,
                "region_class": region_info["region_class"],
                "length_bp": length_bp,
                "mean_depth": mean_depth,
                "breadth": breadth,
                "call": call,
            }
        )
    return rows


def write_pav_calls(rows: list[dict[str, object]], output: Path) -> None:
    columns = [
        "sample_id",
        "region_id",
        "region_class",
        "length_bp",
        "mean_depth",
        "breadth",
        "call",
    ]
    with output.open("w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(columns)
        for row in rows:
            values = [
                f"{row[col]:.4f}" if isinstance(row[col], float) else str(row[col])
                for col in columns
            ]
            writer.writerow(values)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--regions", required=True, type=Path, help="candidate_regions.bed")
    parser.add_argument(
        "--mosdepth-regions", required=True, type=Path, help="mosdepth <prefix>.regions.bed.gz"
    )
    parser.add_argument(
        "--mosdepth-thresholds",
        required=True,
        type=Path,
        help="mosdepth <prefix>.thresholds.bed.gz",
    )
    parser.add_argument("--present-breadth", type=float, default=0.85)
    parser.add_argument("--present-depth", type=float, default=5)
    parser.add_argument("--absent-breadth", type=float, default=0.15)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    region_lengths = read_regions(args.regions)
    mean_depths = read_mean_depth(args.mosdepth_regions)
    threshold_bases = read_threshold_bases(args.mosdepth_thresholds)

    rows = call_pav(
        args.sample_id,
        region_lengths,
        mean_depths,
        threshold_bases,
        args.present_breadth,
        args.present_depth,
        args.absent_breadth,
    )
    write_pav_calls(rows, args.out)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
