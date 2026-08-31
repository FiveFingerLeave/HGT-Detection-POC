#!/usr/bin/env python3
"""Synchronize GFF sequence IDs with a FASTA header id_map (old_id -> new_id)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def load_id_mapping(mapping_tsv: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with mapping_tsv.open("r") as mapfile:
        header = mapfile.readline()
        if header.strip().split("\t") != ["old_id", "new_id"]:
            raise ValueError(f"{mapping_tsv}: expected header 'old_id\\tnew_id'")

        for line_number, line in enumerate(mapfile, start=2):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"{mapping_tsv}:{line_number}: expected two columns")
            old_id, new_id = fields
            mapping[old_id] = new_id
    return mapping


def normalize_gff_seqids(
    input_gff: Path,
    mapping_tsv: Path,
    output_gff: Path,
) -> None:
    """Rewrite the seqid column of every feature line and every
    ``##sequence-region`` pragma using the given old_id -> new_id mapping.
    Fails if a GFF seqid has no entry in the mapping. All other content is
    copied unchanged."""
    mapping = load_id_mapping(mapping_tsv)

    with input_gff.open("r") as infile, output_gff.open("w") as outfile:
        for line_number, line in enumerate(infile, start=1):
            if line.startswith("##sequence-region"):
                parts = line.rstrip("\n").split(" ")
                if len(parts) != 4:
                    raise ValueError(
                        f"{input_gff}:{line_number}: malformed '##sequence-region' pragma"
                    )
                _, old_id, start, end = parts
                new_id = _lookup(old_id, mapping, input_gff, line_number)
                outfile.write(f"##sequence-region {new_id} {start} {end}\n")
                continue

            if line.startswith("#") or not line.strip():
                outfile.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(
                    f"{input_gff}:{line_number}: expected 9 tab-separated GFF columns, "
                    f"got {len(fields)}"
                )

            old_id = fields[0]
            fields[0] = _lookup(old_id, mapping, input_gff, line_number)
            outfile.write("\t".join(fields) + "\n")


def _lookup(
    old_id: str, mapping: dict[str, str], input_gff: Path, line_number: int
) -> str:
    try:
        return mapping[old_id]
    except KeyError:
        raise ValueError(
            f"{input_gff}:{line_number}: seqid '{old_id}' has no entry in the id mapping"
        ) from None


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff", required=True, type=Path, help="Input GFF file")
    parser.add_argument(
        "--mapping",
        required=True,
        type=Path,
        help="old_id/new_id TSV mapping, as produced by normalize_fasta_headers.py",
    )
    parser.add_argument(
        "--output-gff", required=True, type=Path, help="Path for normalized GFF output"
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_gff.parent.mkdir(parents=True, exist_ok=True)

    try:
        normalize_gff_seqids(args.gff, args.mapping, args.output_gff)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
