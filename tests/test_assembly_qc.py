from pathlib import Path

from summarize_assembly_qc import (
    parse_busco_summary,
    parse_quast_report,
    summarize_assembly_qc,
)


def test_parse_quast_report_transposes_correctly(tmp_path: Path) -> None:
    report = tmp_path / "report.tsv"
    report.write_text(
        "Assembly\tISO1\tISO2\n"
        "# contigs\t7\t9\n"
        "Total length\t42000000\t41000000\n"
        "Largest contig\t8800000\t8300000\n"
        "N50\t6000000\t5800000\n"
        "L50\t3\t3\n"
        "GC (%)\t50.1\t50.9\n"
        "Some other metric\tfoo\tbar\n"
    )

    metrics = parse_quast_report(report)

    assert metrics["ISO1"]["num_contigs"] == "7"
    assert metrics["ISO1"]["total_length_bp"] == "42000000"
    assert metrics["ISO2"]["n50"] == "5800000"
    assert "Some other metric" not in metrics["ISO1"]


def test_parse_busco_summary_extracts_percentages(tmp_path: Path) -> None:
    summary = tmp_path / "short_summary.txt"
    summary.write_text(
        "# BUSCO version is: 5.7.1\n"
        "# Summarized benchmarking in BUSCO notation for file genome.fna\n"
        "\tC:98.5%[S:98.0%,D:0.5%],F:0.3%,M:1.2%,n:3817\n"
        "\t3760\tComplete BUSCOs (C)\n"
    )

    result = parse_busco_summary(summary)

    assert result["busco_complete_pct"] == "98.5"
    assert result["busco_single_pct"] == "98.0"
    assert result["busco_duplicated_pct"] == "0.5"
    assert result["busco_fragmented_pct"] == "0.3"
    assert result["busco_missing_pct"] == "1.2"
    assert result["busco_n_markers"] == "3817"


def test_summarize_combines_quast_and_busco(tmp_path: Path) -> None:
    report = tmp_path / "report.tsv"
    report.write_text(
        "Assembly\tISO1\n"
        "# contigs\t7\n"
        "Total length\t42000000\n"
        "N50\t6000000\n"
    )
    summary = tmp_path / "iso1_summary.txt"
    summary.write_text("\tC:98.5%[S:98.0%,D:0.5%],F:0.3%,M:1.2%,n:3817\n")

    rows = summarize_assembly_qc(
        report,
        busco_summaries={"ISO1": summary},
        tiers={"ISO1": "pilot"},
    )

    assert len(rows) == 1
    assert rows[0]["genome_id"] == "ISO1"
    assert rows[0]["tier"] == "pilot"
    assert rows[0]["num_contigs"] == "7"
    assert rows[0]["busco_complete_pct"] == "98.5"


def test_summarize_handles_genome_missing_busco_summary(tmp_path: Path) -> None:
    report = tmp_path / "report.tsv"
    report.write_text("Assembly\tISO1\n# contigs\t7\n")

    rows = summarize_assembly_qc(report, busco_summaries={}, tiers={})

    assert rows[0]["busco_complete_pct"] == ""
