#!/usr/bin/env python3
"""Pivot the combined PAV candidate table (long format: one row per
sample x region) into wide 0/1 matrices (rows = isolates, columns =
regions), split by `region_class` into a core-marker matrix and a
candidate-region matrix (POC guideline, Section 3).

"present" -> 1, "absent" -> 0, "uncertain" -> NaN (excluded from
clustering rather than guessed, per the guideline's own caution that
"uncertain" is better than a wrong "absent" under assembly fragmentation).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

CALL_TO_VALUE = {"present": 1, "absent": 0, "uncertain": float("nan")}


def build_pav_matrix(candidate_table: pd.DataFrame, region_class: str) -> pd.DataFrame:
    subset = candidate_table[candidate_table["region_class"] == region_class]
    subset = subset.assign(value=subset["call"].map(CALL_TO_VALUE))
    return subset.pivot(index="sample_id", columns="region_id", values="value")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-table", required=True, type=Path)
    parser.add_argument("--output-core", required=True, type=Path)
    parser.add_argument("--output-candidate", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_core.parent.mkdir(parents=True, exist_ok=True)
    args.output_candidate.parent.mkdir(parents=True, exist_ok=True)

    table = pd.read_csv(args.candidate_table, sep="\t")

    build_pav_matrix(table, "core").to_csv(args.output_core, sep="\t")
    build_pav_matrix(table, "candidate").to_csv(args.output_candidate, sep="\t")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
