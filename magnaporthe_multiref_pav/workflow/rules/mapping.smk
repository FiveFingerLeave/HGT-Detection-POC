# Phase V (Section 10): Long-read mapping of the pilot test isolates
# against the analytical panel (minimap2), incl. secondary alignments
# for homology evidence and MAPQ-filtered alignments for quantitative
# PAV (Section 10.2, implemented in pav.smk).
#
# Preset choice: The document distinguishes only ONT (map-ont) and PacBio
# HiFi (map-hifi). However, our actual PacBio isolates (B71, K23_123,
# E34) are older PacBio RS II/Sequel RAW subreads (CLR), not
# HiFi/CCS reads - for these, minimap2's own "map-pb" preset is correct,
# not "map-hifi" (see config/samples.tsv: minimap2_preset column,
# set per isolate based on the actual SRA platform).
#
# Only the 5 of 10 pilot isolates with actually downloaded
# raw reads (mappable_sample_ids, Snakefile) are processed here.


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
        # -m 512M explicitly caps samtools sort's memory PER thread -
        # without this, `sort -@ N` alone can already reserve N*768M
        # (default); with several concurrent minimap2_map jobs
        # (Snakemake parallelization across samples) this previously
        # triggered an OOM kill (see docs/decisions.md), even though
        # minimap2 itself only needed ~700 MB per process.
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
