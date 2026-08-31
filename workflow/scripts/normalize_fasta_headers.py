#!/usr/bin/env python3
"""Normalize FASTA headers to <isolate_id><separator><original_id> and record an ID mapping."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def normalize_fasta_headers(
    input_fasta: Path,
    isolate_id: str,
    output_fasta: Path,
    mapping_tsv: Path,
    separator: str = "_",
) -> None:
    """Rewrite FASTA headers as ``<isolate_id><separator><original_id>`` and
    write an ``old_id -> new_id`` TSV mapping. Sequence lines are copied
    unchanged.

    ``separator`` must not itself occur inside ``isolate_id`` (e.g. NCBI
    accessions like ``GCA_004346965.1`` already contain '_'), otherwise
    downstream tools that split on the separator to recover the isolate ID
    cannot do so unambiguously.
    """
    if separator in isolate_id:
        raise ValueError(
            f"separator '{separator}' occurs inside isolate_id '{isolate_id}'; "
            "choose a separator that is not part of the isolate ID"
        )

    seen_new_ids: set[str] = set()
    mapping: list[tuple[str, str]] = []

    with input_fasta.open("r") as infile, output_fasta.open("w") as outfile:
        for line_number, line in enumerate(infile, start=1):
            if not line.startswith(">"):
                outfile.write(line)
                continue

            header = line[1:].rstrip("\n")
            old_id = header.split(maxsplit=1)[0] if header.strip() else ""

            if not old_id:
                raise ValueError(f"{input_fasta}:{line_number}: empty FASTA header")

            new_id = f"{isolate_id}{separator}{old_id}"

            if new_id in seen_new_ids:
                raise ValueError(
                    f"{input_fasta}:{line_number}: duplicate header '{old_id}' "
                    f"produces duplicate ID '{new_id}'"
                )
            seen_new_ids.add(new_id)

            mapping.append((old_id, new_id))
            outfile.write(f">{new_id}\n")

    with mapping_tsv.open("w") as mapfile:
        mapfile.write("old_id\tnew_id\n")
        for old_id, new_id in mapping:
            mapfile.write(f"{old_id}\t{new_id}\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path, help="Input FASTA file")
    parser.add_argument(
        "--isolate-id", required=True, help="Isolate identifier, e.g. GCA_004346965.1"
    )
    parser.add_argument(
        "--output-fasta", required=True, type=Path, help="Path for normalized FASTA output"
    )
    parser.add_argument(
        "--output-mapping",
        required=True,
        type=Path,
        help="Path for old_id/new_id TSV mapping output",
    )
    parser.add_argument(
        "--separator",
        default="_",
        help="Character(s) separating isolate_id from the original id (default: '_')",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.output_mapping.parent.mkdir(parents=True, exist_ok=True)

    try:
        normalize_fasta_headers(
            args.fasta,
            args.isolate_id,
            args.output_fasta,
            args.output_mapping,
            args.separator,
        )
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
