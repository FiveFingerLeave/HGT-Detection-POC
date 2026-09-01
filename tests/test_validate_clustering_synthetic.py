from validate_clustering_synthetic import generate_synthetic_matrices, validate


def test_generate_synthetic_matrices_shapes_and_reproducibility():
    core1, cand1, labels1, recipients1 = generate_synthetic_matrices(
        n_lineages=3, n_isolates_per_lineage=5, n_core_markers=20, n_candidate_regions=4, seed=1
    )
    core2, cand2, labels2, recipients2 = generate_synthetic_matrices(
        n_lineages=3, n_isolates_per_lineage=5, n_core_markers=20, n_candidate_regions=4, seed=1
    )

    assert core1.shape == (15, 20)
    assert cand1.shape == (15, 4)
    assert labels1.equals(labels2)
    assert core1.equals(core2)
    assert recipients1 == recipients2  # same seed -> reproducible


def test_core_matrix_alone_gives_near_perfect_ari():
    result = validate(n_lineages=3, n_isolates_per_lineage=8, n_repeats=5, seed=0)

    assert result["core_ari_mean"] > 0.95
    assert result["n_isolates"] == 24
