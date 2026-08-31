from pathlib import Path

from define_mchr_candidates import classify_contigs, read_rows, write_rows

DEFAULTS = dict(
    core_min_length_bp=1_000_000,
    small_max_length_bp=500_000,
    gc_deviation_threshold=0.03,
    repeat_enrichment_threshold=0.10,
    gene_density_ratio_threshold=0.5,
)


def make_row(**overrides):
    row = {
        "isolate_id": "ISO1",
        "contig_id": "c1",
        "length_bp": 5_000_000,
        "gc_fraction": 0.50,
        "gene_count": 1000,
        "repeat_fraction": 0.10,
        "assembly_unit": "Primary Assembly",
    }
    row.update(overrides)
    return row


def test_large_contig_is_core_like():
    rows = [make_row(contig_id="c1", length_bp=5_000_000)]
    result = classify_contigs(rows, **DEFAULTS)
    assert result[0]["classification"] == "core_like"


def test_non_nuclear_contig_is_excluded_regardless_of_size():
    rows = [make_row(contig_id="mt", length_bp=35_000, assembly_unit="non-nuclear")]
    result = classify_contigs(rows, **DEFAULTS)
    assert result[0]["classification"] == "excluded_non_nuclear"


def test_gray_zone_length_is_uncertain():
    rows = [
        make_row(contig_id="core", length_bp=5_000_000),
        make_row(contig_id="gray", length_bp=700_000, gc_fraction=0.50, repeat_fraction=0.10),
    ]
    result = classify_contigs(rows, **DEFAULTS)
    gray = next(r for r in result if r["contig_id"] == "gray")
    assert gray["classification"] == "uncertain"


def test_small_contig_with_two_evidence_lines_is_mchr_candidate():
    rows = [
        make_row(contig_id="core", length_bp=5_000_000, gc_fraction=0.50, repeat_fraction=0.10),
        make_row(
            contig_id="small",
            length_bp=200_000,
            gc_fraction=0.40,  # deviates by 0.10, above 0.03 threshold
            repeat_fraction=0.30,  # +0.20 over core, above 0.10 threshold
            gene_count=1000,  # same density as core -> not gene-poor
        ),
    ]
    result = classify_contigs(rows, **DEFAULTS)
    small = next(r for r in result if r["contig_id"] == "small")
    assert small["classification"] == "mChr_candidate"


def test_small_contig_with_one_evidence_line_is_accessory_region():
    rows = [
        make_row(contig_id="core", length_bp=5_000_000, gc_fraction=0.50, repeat_fraction=0.10),
        make_row(
            contig_id="small",
            length_bp=200_000,
            gc_fraction=0.50,  # matches core
            repeat_fraction=0.30,  # only this line deviates
            gene_count=1000,
        ),
    ]
    result = classify_contigs(rows, **DEFAULTS)
    small = next(r for r in result if r["contig_id"] == "small")
    assert small["classification"] == "accessory_region"


def test_small_contig_with_two_evidence_lines_and_both_telomeres_is_high_confidence():
    rows = [
        make_row(contig_id="core", length_bp=5_000_000, gc_fraction=0.50, repeat_fraction=0.10),
        make_row(
            contig_id="small",
            length_bp=200_000,
            gc_fraction=0.40,
            repeat_fraction=0.30,
            gene_count=1000,
            telomere_start=True,
            telomere_end=True,
        ),
    ]
    result = classify_contigs(rows, **DEFAULTS)
    small = next(r for r in result if r["contig_id"] == "small")
    assert small["classification"] == "high_confidence_mChr"


def test_small_contig_with_two_evidence_lines_and_one_telomere_stays_candidate():
    rows = [
        make_row(contig_id="core", length_bp=5_000_000, gc_fraction=0.50, repeat_fraction=0.10),
        make_row(
            contig_id="small",
            length_bp=200_000,
            gc_fraction=0.40,
            repeat_fraction=0.30,
            gene_count=1000,
            telomere_start=True,
            telomere_end=False,
        ),
    ]
    result = classify_contigs(rows, **DEFAULTS)
    small = next(r for r in result if r["contig_id"] == "small")
    assert small["classification"] == "mChr_candidate"


def test_small_contig_with_no_evidence_is_uncertain():
    rows = [
        make_row(contig_id="core", length_bp=5_000_000, gc_fraction=0.50, repeat_fraction=0.10),
        make_row(contig_id="small", length_bp=200_000, gc_fraction=0.50, repeat_fraction=0.10),
    ]
    result = classify_contigs(rows, **DEFAULTS)
    small = next(r for r in result if r["contig_id"] == "small")
    assert small["classification"] == "uncertain"


def test_no_core_reference_available_is_uncertain():
    rows = [make_row(contig_id="only_small", length_bp=200_000)]
    result = classify_contigs(rows, **DEFAULTS)
    assert result[0]["classification"] == "uncertain"


def test_missing_gene_count_does_not_crash_and_skips_that_evidence_line():
    rows = [
        make_row(contig_id="core", length_bp=5_000_000, gene_count=None),
        make_row(
            contig_id="small",
            length_bp=200_000,
            gc_fraction=0.40,
            repeat_fraction=0.30,
            gene_count=None,
        ),
    ]
    result = classify_contigs(rows, **DEFAULTS)
    small = next(r for r in result if r["contig_id"] == "small")
    assert small["classification"] == "mChr_candidate"


def test_read_write_round_trip(tmp_path: Path) -> None:
    input_tsv = tmp_path / "in.tsv"
    input_tsv.write_text(
        "isolate_id\tcontig_id\tlength_bp\tgc_fraction\tgene_count\trepeat_fraction\tassembly_unit\trole\ttelomere_start\ttelomere_end\n"
        "ISO1\tc1\t5000000\t0.5000\t1000\t0.1000\tPrimary Assembly\tassembled-molecule\tTrue\tTrue\n"
        "ISO1\tmt\t35000\t0.2900\t\t0.0500\tnon-nuclear\tassembled-molecule\tFalse\tFalse\n"
    )

    rows = read_rows(input_tsv)
    assert rows[1]["gene_count"] is None
    assert rows[1]["assembly_unit"] == "non-nuclear"

    classified = classify_contigs(rows, **DEFAULTS)
    output_tsv = tmp_path / "out.tsv"
    write_rows(classified, output_tsv)

    written = output_tsv.read_text().splitlines()
    assert written[0].split("\t")[-1] == "classification"
    assert written[2].endswith("excluded_non_nuclear")
