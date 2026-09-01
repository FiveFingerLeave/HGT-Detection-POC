from validate_clustering_synthetic import generate_synthetic_matrices, validate


def test_generate_synthetic_matrices_shapes_and_reproducibility():
    core1, cand1, labels1, recipients1, hgt_col1 = generate_synthetic_matrices(
        n_lineages=3, n_isolates_per_lineage=5, n_core_markers=20, n_candidate_regions=4, seed=1
    )
    core2, cand2, labels2, recipients2, hgt_col2 = generate_synthetic_matrices(
        n_lineages=3, n_isolates_per_lineage=5, n_core_markers=20, n_candidate_regions=4, seed=1
    )

    assert core1.shape == (15, 20)
    assert cand1.shape == (15, 4)
    assert labels1.equals(labels2)
    assert core1.equals(core2)
    assert recipients1 == recipients2  # same seed -> reproducible
    assert hgt_col1 == hgt_col2 == 0


def test_core_matrix_alone_gives_near_perfect_ari():
    result = validate(n_lineages=3, n_isolates_per_lineage=8, n_repeats=5, seed=0)

    assert result["core_ari_mean"] > 0.95
    assert result["n_isolates"] == 24


def test_per_region_method_has_higher_detection_rate_than_aggregate_method():
    # This directly compares Method A (aggregate whole-matrix ARI) vs.
    # Method B (per-region ARI) on the same synthetic scenario - the
    # motivation for building Method B in the first place (see
    # docs/decisions.md). Note: with 4 lineages but a binary (2-valued)
    # per-region readout, baseline regions cannot reach ARI near 1.0
    # (pigeonhole: at least two lineages always share a bit value by
    # chance) - the meaningful comparison is RELATIVE, injected region vs.
    # baseline regions, not an absolute "near 1.0" expectation.
    result = validate(n_repeats=15, seed=2)

    assert result["per_region_detection_rate"] >= result["fraction_trials_with_clear_drop"]
    assert result["hgt_region_ari_mean"] < result["baseline_region_ari_mean"]
