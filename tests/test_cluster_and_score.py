import numpy as np
import pandas as pd

from cluster_and_score import cluster_and_score, per_region_ari


def test_cluster_and_score_recovers_perfect_lineage_structure():
    # Two clearly separated lineage blocks, each internally identical.
    rng = np.random.default_rng(0)
    lineage_a = pd.DataFrame(np.tile([1, 1, 0, 0], (10, 1)))
    lineage_b = pd.DataFrame(np.tile([0, 0, 1, 1], (10, 1)))
    matrix = pd.concat([lineage_a, lineage_b], ignore_index=True)
    labels = pd.Series(["A"] * 10 + ["B"] * 10)

    _, clusters, ari = cluster_and_score(matrix, labels, k=2, n_components=2)

    assert ari > 0.9


def test_cluster_and_score_handles_too_few_isolates():
    matrix = pd.DataFrame([[1, 0, 1]])
    labels = pd.Series(["A"])

    coords, clusters, ari = cluster_and_score(matrix, labels, k=2)

    assert coords is None
    assert clusters is None
    assert ari != ari  # nan


def test_per_region_ari_high_for_lineage_clean_region():
    # region perfectly separates two lineages -> its own partition IS the
    # lineage partition -> ARI should be 1.0
    labels = pd.Series(["A"] * 5 + ["B"] * 5)
    matrix = pd.DataFrame({"regionX": [1] * 5 + [0] * 5})

    scores = per_region_ari(matrix, labels)

    assert scores["regionX"] == 1.0


def test_per_region_ari_low_for_cross_lineage_sharing():
    # "present" cuts across both lineages instead of respecting them.
    labels = pd.Series(["A"] * 5 + ["B"] * 5)
    matrix = pd.DataFrame({"regionX": [1, 1, 0, 0, 0, 1, 1, 0, 0, 0]})

    scores = per_region_ari(matrix, labels)

    assert scores["regionX"] < 0.3


def test_per_region_ari_nan_for_constant_region():
    labels = pd.Series(["A"] * 5 + ["B"] * 5)
    matrix = pd.DataFrame({"regionX": [1] * 10})

    scores = per_region_ari(matrix, labels)

    assert scores["regionX"] != scores["regionX"]  # nan


def test_per_region_ari_drops_missing_calls():
    labels = pd.Series(["A"] * 5 + ["B"] * 5)
    # region perfectly clean except for one NaN ("uncertain" PAV call),
    # which must be excluded rather than breaking/skewing the score.
    matrix = pd.DataFrame({"regionX": [1, 1, 1, 1, float("nan"), 0, 0, 0, 0, 0]})

    scores = per_region_ari(matrix, labels)

    assert scores["regionX"] == 1.0


def test_per_region_ari_returns_one_score_per_column():
    labels = pd.Series(["A", "A", "B", "B"])
    matrix = pd.DataFrame({"clean": [1, 1, 0, 0], "mixed": [1, 0, 1, 0]})

    scores = per_region_ari(matrix, labels)

    assert list(scores.index) == ["clean", "mixed"]
    assert scores["clean"] > scores["mixed"]
