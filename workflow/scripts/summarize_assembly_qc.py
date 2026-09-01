#!/usr/bin/env python3
"""Combine QUAST contiguity metrics and BUSCO completeness scores into one
per-genome assembly QC summary table.

Establishes a standardized, comparable quality baseline across the pilot
isolates and the reference panel (config/references.yaml) - assembly-level
metrics only, since no raw reads are available for these genomes (they are
already-assembled NCBI submissions, not something this project sequenced
and assembled itself).
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

QUAST_METRICS = {
    "# contigs": "num_contigs",
    "Total length": "total_length_bp",
    "Largest contig": "largest_contig_bp",
    "N50": "n50",
    "L50": "l50",
    "GC (%)": "gc_percent",
}

BUSCO_SUMMARY_PATTERN = re.compile(
    r"C:(?P<complete>[\d.]+)%\[S:(?P<single>[\d.]+)%,D:(?P<duplicated>[\d.]+)%\],"
    r"F:(?P<fragmented>[\d.]+)%,M:(?P<missing>[\d.]+)%,n:(?P<n>\d+)"
)


def parse_quast_report(path: Path) -> dict[str, dict[str, str]]:
    """Parse QUAST's transposed report.tsv: first column is the metric
    name, remaining columns are one per genome (the --labels used)."""
    with path.open("r", newline="") as infile:
        rows = list(csv.reader(infile, delimiter="\t"))

    if not rows:
        return {}

    genome_ids = rows[0][1:]
    metrics: dict[str, dict[str, str]] = {gid: {} for gid in genome_ids}

    for row in rows[1:]:
        if not row:
            continue
        metric_name = row[0]
        if metric_name not in QUAST_METRICS:
            continue
        key = QUAST_METRICS[metric_name]
        for gid, value in zip(genome_ids, row[1:]):
            metrics[gid][key] = value

    return metrics


def parse_busco_summary(path: Path) -> dict[str, str]:
    """Parse a BUSCO short_summary file's 'C:X%[S:Y%,D:Z%],F:A%,M:B%,n:N' line."""
    text = path.read_text()
    match = BUSCO_SUMMARY_PATTERN.search(text)
    if not match:
        raise ValueError(f"{path}: could not find BUSCO summary line")

    return {
        "busco_complete_pct": match.group("complete"),
        "busco_single_pct": match.group("single"),
        "busco_duplicated_pct": match.group("duplicated"),
        "busco_fragmented_pct": match.group("fragmented"),
        "busco_missing_pct": match.group("missing"),
        "busco_n_markers": match.group("n"),
    }


def summarize_assembly_qc(
    quast_report: Path,
    busco_summaries: dict[str, Path],
    tiers: dict[str, str],
) -> list[dict[str, object]]:
    quast_metrics = parse_quast_report(quast_report)

    genome_ids = sorted(set(quast_metrics) | set(busco_summaries))
    rows = []
    for genome_id in genome_ids:
        row: dict[str, object] = {
            "genome_id": genome_id,
            "tier": tiers.get(genome_id, ""),
        }
        row.update({key: "" for key in QUAST_METRICS.values()})
        row.update(quast_metrics.get(genome_id, {}))
        row.update(
            {
                "busco_complete_pct": "",
                "busco_single_pct": "",
                "busco_duplicated_pct": "",
                "busco_fragmented_pct": "",
                "busco_missing_pct": "",
                "busco_n_markers": "",
            }
        )
        if genome_id in busco_summaries:
            row.update(parse_busco_summary(busco_summaries[genome_id]))
        rows.append(row)

    return rows


def write_summary(rows: list[dict[str, object]], output: Path) -> None:
    columns = [
        "genome_id",
        "tier",
        "num_contigs",
        "total_length_bp",
        "largest_contig_bp",
        "n50",
        "l50",
        "gc_percent",
        "busco_complete_pct",
        "busco_single_pct",
        "busco_duplicated_pct",
        "busco_fragmented_pct",
        "busco_missing_pct",
        "busco_n_markers",
    ]
    with output.open("w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([row.get(col, "") for col in columns])


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quast-report", required=True, type=Path)
    parser.add_argument(
        "--busco-summary",
        action="append",
        default=[],
        metavar="GENOME_ID=PATH",
        dest="busco_summary_pairs",
        help="Repeatable: genome_id=path/to/short_summary.txt",
    )
    parser.add_argument(
        "--tier",
        action="append",
        default=[],
        metavar="GENOME_ID=TIER",
        dest="tier_pairs",
        help="Repeatable: genome_id=pilot|structural|supplementary",
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    busco_summaries = {}
    for pair in args.busco_summary_pairs:
        genome_id, path = pair.split("=", 1)
        busco_summaries[genome_id] = Path(path)

    tiers = dict(pair.split("=", 1) for pair in args.tier_pairs)

    rows = summarize_assembly_qc(args.quast_report, busco_summaries, tiers)
    write_summary(rows, args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
