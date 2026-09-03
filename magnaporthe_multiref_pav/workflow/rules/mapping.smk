# Phase V (Section 10): Long-Read-Mapping der Pilot-Testisolate gegen
# das analytische Panel (minimap2), inkl. secondary Alignments fuer
# Homologie-Nachweis und MAPQ-gefilterte Alignments fuer quantitative
# PAV (Abschnitt 10.2, umgesetzt in pav.smk).
#
# Preset-Wahl: Das Dokument unterscheidet nur ONT (map-ont) und PacBio
# HiFi (map-hifi). Unsere tatsaechlichen PacBio-Isolate (B71, K23_123,
# E34) sind aber aeltere PacBio RS II/Sequel RAW-Subreads (CLR), keine
# HiFi/CCS-Reads - dafuer ist minimap2s eigenes "map-pb"-Preset korrekt,
# nicht "map-hifi" (siehe config/samples.tsv: minimap2_preset-Spalte,
# pro Isolat anhand der tatsaechlichen SRA-Plattform gesetzt).
#
# Nur die 5 von 10 Pilotisolaten mit tatsaechlich heruntergeladenen
# Rohreads (mappable_sample_ids, Snakefile) werden hier prozessiert.


rule minimap2_map:
    threads: config["threads_default"]
    input:
        panel="results/panel/Mo_multiref_panel_v1.fa",
        panel_fai="results/panel/Mo_multiref_panel_v1.fa.fai",
        fastq=lambda wc: sample_fastq[wc.sample_id],
    output:
        bam="results/mapping/{sample_id}.panel.bam",
        bai="results/mapping/{sample_id}.panel.bam.bai",
    params:
        preset=lambda wc: sample_preset[wc.sample_id],
    log:
        "logs/mapping/{sample_id}_minimap2.log",
    conda:
        "../../envs/core.yaml"
    shell:
        # -m 512M begrenzt samtools sorts Speicher PRO Thread explizit -
        # ohne das kann `sort -@ N` allein schon N*768M (Default)
        # reservieren; bei mehreren gleichzeitigen minimap2_map-Jobs
        # (Snakemake-Parallelisierung ueber Samples) hat das zuvor einen
        # OOM-Kill ausgeloest (siehe docs/decisions.md), obwohl minimap2
        # selbst pro Prozess nur ~700 MB brauchte.
        "mkdir -p results/mapping logs/mapping && "
        "minimap2 -ax {params.preset} --secondary=yes -t {threads} "
        "{input.panel} {input.fastq} 2> {log} "
        "| samtools sort -@ 2 -m 512M -o {output.bam} - "
        ">> {log} 2>&1 && "
        "samtools index {output.bam}"


rule mapping_flagstat:
    input:
        "results/mapping/{sample_id}.panel.bam",
    output:
        "results/mapping/{sample_id}.flagstat.txt",
    conda:
        "../../envs/core.yaml"
    shell:
        "samtools flagstat {input} > {output}"


rule mapping_coverage_by_contig:
    input:
        "results/mapping/{sample_id}.panel.bam",
    output:
        "results/mapping/{sample_id}.coverage_by_contig.tsv",
    conda:
        "../../envs/core.yaml"
    shell:
        "samtools coverage {input} > {output}"


rule mapping_all:
    input:
        expand("results/mapping/{sample_id}.flagstat.txt", sample_id=mappable_sample_ids),
        expand("results/mapping/{sample_id}.coverage_by_contig.tsv", sample_id=mappable_sample_ids),
