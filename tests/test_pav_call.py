from pathlib import Path

from pav_call import call_pav, read_regions, read_threshold_bases


def test_read_regions_returns_lengths_and_default_class(tmp_path: Path) -> None:
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t1000\t2000\tregionA\nchr1\t5000\t5500\tregionB\n")

    regions = read_regions(bed)

    assert regions == {
        "regionA": {"length_bp": 1000, "region_class": ""},
        "regionB": {"length_bp": 500, "region_class": ""},
    }


def test_read_regions_captures_region_class(tmp_path: Path) -> None:
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t0\t1000\tregionA\tcore\nchr1\t2000\t3000\tregionB\tcandidate\n")

    regions = read_regions(bed)

    assert regions["regionA"]["region_class"] == "core"
    assert regions["regionB"]["region_class"] == "candidate"


def test_call_pav_present_when_breadth_and_depth_pass():
    rows = call_pav(
        "sample1",
        regions={"regionA": {"length_bp": 1000, "region_class": "candidate"}},
        mean_depths={"regionA": 10.0},
        threshold_bases={"regionA": 900},  # 90% breadth at present_depth
        present_breadth=0.85,
        present_depth=5,
        absent_breadth=0.15,
    )

    assert rows[0]["call"] == "present"
    assert rows[0]["breadth"] == 0.9
    assert rows[0]["region_class"] == "candidate"


def test_call_pav_absent_when_breadth_low():
    rows = call_pav(
        "sample1",
        regions={"regionA": {"length_bp": 1000, "region_class": "candidate"}},
        mean_depths={"regionA": 0.5},
        threshold_bases={"regionA": 50},  # 5% breadth
        present_breadth=0.85,
        present_depth=5,
        absent_breadth=0.15,
    )

    assert rows[0]["call"] == "absent"


def test_call_pav_uncertain_in_gray_zone():
    rows = call_pav(
        "sample1",
        regions={"regionA": {"length_bp": 1000, "region_class": "candidate"}},
        mean_depths={"regionA": 3.0},
        threshold_bases={"regionA": 500},  # 50% breadth: neither present nor absent
        present_breadth=0.85,
        present_depth=5,
        absent_breadth=0.15,
    )

    assert rows[0]["call"] == "uncertain"


def test_call_pav_uncertain_when_breadth_high_but_depth_low():
    rows = call_pav(
        "sample1",
        regions={"regionA": {"length_bp": 1000, "region_class": "candidate"}},
        mean_depths={"regionA": 2.0},  # below present_depth
        threshold_bases={"regionA": 950},  # 95% breadth
        present_breadth=0.85,
        present_depth=5,
        absent_breadth=0.15,
    )

    assert rows[0]["call"] == "uncertain"


def test_call_pav_missing_region_data_defaults_to_absent():
    # A region with no coverage data at all (e.g. no reads mapped there)
    rows = call_pav(
        "sample1",
        regions={"regionA": {"length_bp": 1000, "region_class": "candidate"}},
        mean_depths={},
        threshold_bases={},
        present_breadth=0.85,
        present_depth=5,
        absent_breadth=0.15,
    )

    assert rows[0]["call"] == "absent"
    assert rows[0]["breadth"] == 0.0


def test_read_threshold_bases_skips_header(tmp_path: Path) -> None:
    bed = tmp_path / "thresholds.bed"
    bed.write_text("#chrom\tstart\tend\tregion\t5X\nchr1\t0\t1000\tregionA\t900\n")

    bases = read_threshold_bases(bed)

    assert bases == {"regionA": 900}
