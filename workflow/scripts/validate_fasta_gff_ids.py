#!/usr/bin/env python3
"""Validate that every GFF seqid has a matching FASTA sequence."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def read_fasta_ids(fasta: Path) -> set[str]:
    ids: set[str] = set()
    with fasta.open("r") as infile:
        for line in infile:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                ids.add(header.split(maxsplit=1)[0] if header.strip() else "")
    return ids


def read_gff_seqids(gff: Path) -> set[str]:
    ids: set[str] = set()
    with gff.open("r") as infile:
        for line in infile:
            if line.startswith("##sequence-region"):
                parts = line.rstrip("\n").split(" ")
                if len(parts) == 4:
                    ids.add(parts[1])
                continue

            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if fields:
                ids.add(fields[0])
    return ids


def validate_fasta_gff_ids(fasta: Path, gff: Path) -> list[str]:
    """Return the sorted list of GFF seqids that have no matching FASTA
    sequence. An empty list means the files are compatible."""
    fasta_ids = read_fasta_ids(fasta)
    gff_ids = read_gff_seqids(gff)
    return sorted(gff_ids - fasta_ids)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path, help="FASTA file")
    parser.add_argument("--gff", required=True, type=Path, help="GFF file")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    missing = validate_fasta_gff_ids(args.fasta, args.gff)

    if missing:
        print(
            f"Error: {len(missing)} GFF seqid(s) have no matching FASTA sequence:",
            file=sys.stderr,
        )
        for seqid in missing:
            print(f"  {seqid}", file=sys.stderr)
        return 1

    print("OK: all GFF seqids have a matching FASTA sequence.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
