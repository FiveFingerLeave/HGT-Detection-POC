#!/usr/bin/env python3
"""PCA + k-means clustering of a PAV matrix, scored against known host-
lineage labels via Adjusted Rand Index (ARI). POC guideline Section 3:
core-genome markers should cluster isolates by host lineage (high ARI);
HGT-candidate regions should show discordant, lineage-crossing clustering
(lower ARI) - the central evidence signature for horizontal transfer.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score


def per_region_ari(pav_matrix: pd.DataFrame, host_labels: pd.Series) -> pd.Series:
    """ARI between each single region's presence/absence and host lineage,
    independent of every other region - a targeted, per-locus alternative
    to `cluster_and_score`'s whole-matrix PCA/k-means ARI.

    A region that is cleanly inherited within lineages (present in some
    lineages, absent in others, but uniform WITHIN each lineage) gets a
    high ARI - its own presence/absence already IS a partition that lines
    up with host lineage. A region shared across lineages in a way that
    doesn't respect lineage boundaries (the HGT signature) gets a low ARI.
    Isolates with a missing/NaN call for a given region are dropped from
    that region's score only. A region with no variance (present or
    absent everywhere, once NaNs are dropped) or too few remaining
    isolates has an undefined ARI and is reported as NaN.
    """
    scores: dict[object, float] = {}
    for region_id in pav_matrix.columns:
        column = pav_matrix[region_id]
        valid = column.notna()
        if valid.sum() < 2 or column[valid].nunique() < 2:
            scores[region_id] = float("nan")
            continue
        scores[region_id] = adjusted_rand_score(host_labels[valid], column[valid])
    return pd.Series(scores, name="ari")


def cluster_and_score(
    pav_matrix: pd.DataFrame, host_labels: pd.Series, k: int, n_components: int = 10
) -> tuple["pd.Series | None", "pd.Series | None", float]:
    """Return (pca_coords, cluster_labels, ari). None/nan if the matrix has
    too few isolates or regions to cluster meaningfully."""
    n_isolates, n_regions = pav_matrix.shape
    if n_isolates < 2 or n_regions < 1:
        return None, None, float("nan")

    max_components = min(n_components, n_regions, n_isolates - 1)
    if max_components < 1:
        return None, None, float("nan")

    pca = PCA(n_components=max_components)
    coords = pca.fit_transform(pav_matrix.values)

    k = min(k, n_isolates)
    clusters = KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(coords)
    ari = adjusted_rand_score(host_labels, clusters)

    return coords, clusters, ari


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core", required=True, type=Path, help="Core PAV matrix TSV")
    parser.add_argument(
        "--candidate", required=True, type=Path, help="Candidate-region PAV matrix TSV"
    )
    parser.add_argument(
        "--host-labels",
        required=True,
        type=Path,
        help="TSV with sample_id and host_lineage columns (e.g. config/samples.tsv)",
    )
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument("--n-components", type=int, default=10)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--per-region-output",
        type=Path,
        default=None,
        help="Optional: write per-region ARI (candidate matrix only) to this TSV - "
        "a more sensitive, targeted complement to the whole-matrix ARI above "
        "(see docs/decisions.md).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    labels_df = pd.read_csv(args.host_labels, sep="\t", dtype=str).set_index("sample_id")

    rows = []
    candidate_matrix = None
    candidate_host_labels = None
    for name, path in (("core", args.core), ("candidate", args.candidate)):
        matrix = pd.read_csv(path, sep="\t", index_col=0)
        host_labels = labels_df.loc[matrix.index, "host_lineage"]
        _, _, ari = cluster_and_score(matrix, host_labels, args.k, args.n_components)
        rows.append(
            {
                "matrix": name,
                "n_isolates": matrix.shape[0],
                "n_regions": matrix.shape[1],
                "k": args.k,
                "ari": ari,
            }
        )
        if name == "candidate":
            candidate_matrix, candidate_host_labels = matrix, host_labels

    with args.output.open("w", newline="") as outfile:
        writer = csv.DictWriter(
            outfile, fieldnames=["matrix", "n_isolates", "n_regions", "k", "ari"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    if args.per_region_output is not None:
        args.per_region_output.parent.mkdir(parents=True, exist_ok=True)
        per_region = per_region_ari(candidate_matrix, candidate_host_labels)
        per_region.rename_axis("region_id").reset_index(name="ari").to_csv(
            args.per_region_output, sep="\t", index=False
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
