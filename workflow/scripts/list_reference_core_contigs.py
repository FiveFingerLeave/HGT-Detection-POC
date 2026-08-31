#!/usr/bin/env python3
"""List the "core" contig IDs of a reference genome: nuclear contigs at/above
--core-min-length-bp, minus any explicitly excluded (documented accessory /
mini-chromosome) contigs.

Length, not the NCBI sequence report's `role` field, is used as the primary
signal: contig-level (un-scaffolded) reference assemblies report every
contig as 'unplaced-scaffold' regardless of size, so `role` alone cannot
distinguish a multi-Mb chromosome from a small fragment for those
assemblies. `assemblyUnit` (non-nuclear / mitochondrial) is still used when
available, since it is reliable across assembly types.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def read_contig_lengths(fasta: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    contig_id: str | None = None
    with fasta.open("r") as infile:
        for line in infile:
            if line.startswith(">"):
                contig_id = line[1:].split(maxsplit=1)[0].rstrip("\n")
                lengths[contig_id] = 0
                continue
            if contig_id is not None:
                lengths[contig_id] += len(line.strip())
    return lengths


def read_non_nuclear_ids(sequence_report: Path | None) -> set[str]:
    if sequence_report is None:
        return set()
    non_nuclear = set()
    with sequence_report.open("r") as infile:
        for line in infile:
            if not line.strip():
                continue
            record = json.loads(line)
            if record.get("assemblyUnit") == "non-nuclear":
                non_nuclear.add(record.get("genbankAccession"))
    return non_nuclear


def list_core_contigs(
    fasta: Path,
    core_min_length_bp: int,
    sequence_report: Path | None = None,
    exclude: set[str] | None = None,
) -> list[str]:
    lengths = read_contig_lengths(fasta)
    non_nuclear = read_non_nuclear_ids(sequence_report)
    exclude = exclude or set()

    return [
        contig_id
        for contig_id, length in lengths.items()
        if length >= core_min_length_bp
        and contig_id not in non_nuclear
        and contig_id not in exclude
    ]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--core-min-length-bp", type=int, default=4_000_000)
    parser.add_argument("--sequence-report", type=Path, default=None)
    parser.add_argument(
        "--exclude",
        default="",
        help="Comma-separated contig IDs to exclude (documented accessory contigs)",
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    exclude = {c.strip() for c in args.exclude.split(",") if c.strip()}
    core_contigs = list_core_contigs(
        args.fasta, args.core_min_length_bp, args.sequence_report, exclude
    )

    with args.output.open("w") as outfile:
        for contig_id in core_contigs:
            outfile.write(f"{contig_id}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
