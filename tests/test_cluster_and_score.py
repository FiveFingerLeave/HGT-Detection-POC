import numpy as np
import pandas as pd

from cluster_and_score import cluster_and_score


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
