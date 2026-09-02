#!/usr/bin/env python3
"""Rename FASTA contig headers to '>{genome_id}__{original_contig_name}'
(Dokumentation/multireferenzpanel_pav_workflow.md, Section 5.3) so contig
names are globally unique across the whole reference panel. Writes a
name-mapping table alongside the renamed copy; the original FASTA is
never modified.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def rename_headers(in_fasta: Path, genome_id: str, out_fasta: Path) -> list[tuple[str, str]]:
    mapping = []
    with in_fasta.open() as infile, out_fasta.open("w") as outfile:
        for line in infile:
            if line.startswith(">"):
                original = line[1:].split()[0].strip()
                new_name = f"{genome_id}__{original}"
                mapping.append((genome_id, original, new_name))
                outfile.write(f">{new_name}\n")
            else:
                outfile.write(line)
    return mapping


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--genome-id", required=True)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument(
        "--output-map",
        required=True,
        type=Path,
        help="TSV appended to (created with header if absent): genome_id, original_contig, new_contig",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.output_map.parent.mkdir(parents=True, exist_ok=True)

    mapping = rename_headers(args.fasta, args.genome_id, args.output_fasta)

    write_header = not args.output_map.exists()
    with args.output_map.open("a", newline="") as mapfile:
        writer = csv.writer(mapfile, delimiter="\t")
        if write_header:
            writer.writerow(["genome_id", "original_contig", "new_contig"])
        writer.writerows(mapping)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
