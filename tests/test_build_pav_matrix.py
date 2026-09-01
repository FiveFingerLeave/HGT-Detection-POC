import pandas as pd

from build_pav_matrix import build_pav_matrix


def test_build_pav_matrix_splits_by_region_class_and_maps_calls():
    table = pd.DataFrame(
        [
            {"sample_id": "iso1", "region_id": "core1", "region_class": "core", "call": "present"},
            {"sample_id": "iso1", "region_id": "cand1", "region_class": "candidate", "call": "absent"},
            {"sample_id": "iso2", "region_id": "core1", "region_class": "core", "call": "absent"},
            {"sample_id": "iso2", "region_id": "cand1", "region_class": "candidate", "call": "uncertain"},
        ]
    )

    core = build_pav_matrix(table, "core")
    candidate = build_pav_matrix(table, "candidate")

    assert core.loc["iso1", "core1"] == 1
    assert core.loc["iso2", "core1"] == 0
    assert candidate.loc["iso1", "cand1"] == 0
    assert pd.isna(candidate.loc["iso2", "cand1"])
