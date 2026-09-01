#!/usr/bin/env python3
"""Unit test of the clustering METHOD (not a biological result), per the
POC guideline Section 3 "Empfehlung: Methodenvalidierung an synthetischen
Daten (jetzt schon möglich)". Runs independently of any real sequence
data.

Ground truth: `n_lineages` host lineages, each with its own fixed
core-marker profile (isolates of the same lineage are core-identical,
different lineages differ) and its own fixed candidate-region baseline
profile. Into the candidate matrix, one region is then forced "present" in
a subset of isolates drawn ACROSS lineages (simulated horizontal
transfer), regardless of their lineage's baseline.

Expected result (the guideline's own success criteria):
- Core matrix, no injected HGT -> ARI vs. true lineage labels near 1.
- Candidate matrix, with injected cross-lineage sharing -> ARI
  significantly lower than the core ARI.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from cluster_and_score import cluster_and_score


def generate_synthetic_matrices(
    n_lineages: int = 4,
    n_isolates_per_lineage: int = 10,
    n_core_markers: int = 200,
    n_candidate_regions: int = 30,
    hgt_fraction: float = 0.15,
    seed: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, list[str]]:
    rng = np.random.default_rng(seed)

    lineage_names = [f"lineage_{i}" for i in range(n_lineages)]
    isolate_ids = [
        f"{lineage}_iso{j}"
        for lineage in lineage_names
        for j in range(n_isolates_per_lineage)
    ]
    host_labels = pd.Series(
        [iso.rsplit("_iso", 1)[0] for iso in isolate_ids], index=isolate_ids
    )

    # Core matrix: one fixed binary profile per lineage, shared exactly by
    # all its isolates -> should cluster perfectly by lineage.
    core_profiles = {
        lineage: rng.integers(0, 2, size=n_core_markers) for lineage in lineage_names
    }
    core_matrix = pd.DataFrame(
        [core_profiles[host_labels[iso]] for iso in isolate_ids], index=isolate_ids
    )

    # Candidate matrix: same idea as a baseline...
    candidate_profiles = {
        lineage: rng.integers(0, 2, size=n_candidate_regions) for lineage in lineage_names
    }
    candidate_matrix = pd.DataFrame(
        [candidate_profiles[host_labels[iso]] for iso in isolate_ids], index=isolate_ids
    )

    # ...but inject cross-lineage sharing of ONE candidate region into a
    # subset of isolates drawn from every lineage, simulating HGT.
    hgt_region_col = 0
    n_recipients = max(1, int(round(hgt_fraction * len(isolate_ids))))
    hgt_recipients = list(
        rng.choice(isolate_ids, size=n_recipients, replace=False)
    )
    candidate_matrix.loc[hgt_recipients, hgt_region_col] = 1

    return core_matrix, candidate_matrix, host_labels, hgt_recipients


def _validate_once(
    n_lineages: int,
    n_isolates_per_lineage: int,
    n_core_markers: int,
    n_candidate_regions: int,
    hgt_fraction: float,
    k: int,
    seed: int,
) -> dict[str, float]:
    core_matrix, candidate_matrix, host_labels, hgt_recipients = generate_synthetic_matrices(
        n_lineages, n_isolates_per_lineage, n_core_markers, n_candidate_regions, hgt_fraction, seed
    )

    _, _, core_ari = cluster_and_score(core_matrix, host_labels, k)
    _, cand_clusters, candidate_ari = cluster_and_score(candidate_matrix, host_labels, k)

    cand_cluster_series = pd.Series(cand_clusters, index=candidate_matrix.index)
    recipient_clusters = cand_cluster_series.loc[hgt_recipients]
    largest_shared_cluster_fraction = (
        recipient_clusters.value_counts(normalize=True).iloc[0]
        if len(recipient_clusters)
        else float("nan")
    )

    return {
        "core_ari": core_ari,
        "candidate_ari": candidate_ari,
        "hgt_recipients_grouped_fraction": largest_shared_cluster_fraction,
    }


def validate(
    n_lineages: int = 4,
    n_isolates_per_lineage: int = 10,
    n_core_markers: int = 200,
    n_candidate_regions: int = 5,
    hgt_fraction: float = 0.3,
    k: int | None = None,
    n_repeats: int = 20,
    seed: int = 0,
) -> dict[str, object]:
    """Monte-Carlo validation over `n_repeats` random draws.

    A SINGLE draw is unreliable: with only a handful of candidate regions,
    whether one injected discordant region visibly perturbs the whole-
    matrix PCA/k-means clustering depends on how it happens to align with
    the other (unrelated, randomly-drawn) candidate baselines - empirically
    ~35-65% of single seeds fail to show a clear ARI drop even though the
    injected signal is real (see docs/decisions.md). Averaging over many
    repeats is the honest way to characterize the method's sensitivity,
    rather than reporting one favorable seed.
    """
    k = k or n_lineages
    rng = np.random.default_rng(seed)
    trial_seeds = rng.integers(0, 2**31 - 1, size=n_repeats)

    trials = [
        _validate_once(
            n_lineages,
            n_isolates_per_lineage,
            n_core_markers,
            n_candidate_regions,
            hgt_fraction,
            k,
            int(trial_seed),
        )
        for trial_seed in trial_seeds
    ]
    trials_df = pd.DataFrame(trials)

    n_isolates = n_lineages * n_isolates_per_lineage
    n_hgt_recipients = max(1, int(round(hgt_fraction * n_isolates)))

    return {
        "core_ari_mean": trials_df["core_ari"].mean(),
        "core_ari_std": trials_df["core_ari"].std(),
        "candidate_ari_mean": trials_df["candidate_ari"].mean(),
        "candidate_ari_std": trials_df["candidate_ari"].std(),
        "ari_dropped_mean": trials_df["core_ari"].mean() - trials_df["candidate_ari"].mean(),
        "fraction_trials_with_clear_drop": (
            (trials_df["core_ari"] - trials_df["candidate_ari"]) > 0.2
        ).mean(),
        "hgt_recipients_grouped_fraction_mean": trials_df[
            "hgt_recipients_grouped_fraction"
        ].mean(),
        "n_isolates": n_isolates,
        "n_hgt_recipients": n_hgt_recipients,
        "n_repeats": n_repeats,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-lineages", type=int, default=4)
    parser.add_argument("--n-isolates-per-lineage", type=int, default=10)
    parser.add_argument("--n-core-markers", type=int, default=200)
    parser.add_argument("--n-candidate-regions", type=int, default=5)
    parser.add_argument("--hgt-fraction", type=float, default=0.3)
    parser.add_argument("--n-repeats", type=int, default=20)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    result = validate(
        args.n_lineages,
        args.n_isolates_per_lineage,
        args.n_core_markers,
        args.n_candidate_regions,
        args.hgt_fraction,
        n_repeats=args.n_repeats,
        seed=args.seed,
    )

    print(f"Monte-Carlo validation over {result['n_repeats']} random draws:")
    print(
        f"Core ARI (no injected HGT):        "
        f"{result['core_ari_mean']:.3f} +/- {result['core_ari_std']:.3f}  (expect: near 1.0)"
    )
    print(
        f"Candidate ARI (with injected HGT):  "
        f"{result['candidate_ari_mean']:.3f} +/- {result['candidate_ari_std']:.3f}  (expect: << core ARI)"
    )
    print(f"Mean ARI drop:                       {result['ari_dropped_mean']:.3f}")
    print(
        f"Trials with a clear ARI drop (>0.2): "
        f"{result['fraction_trials_with_clear_drop']:.1%}"
    )
    print(
        f"HGT recipients grouped together:     "
        f"{result['hgt_recipients_grouped_fraction_mean']:.1%} avg. in the same candidate cluster "
        f"({result['n_hgt_recipients']}/{result['n_isolates']} isolates)"
    )

    ok = result["core_ari_mean"] > 0.9 and result["ari_dropped_mean"] > 0.2
    print("\nPASS" if ok else "\nFAIL: clustering method did not reproduce the expected pattern")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
