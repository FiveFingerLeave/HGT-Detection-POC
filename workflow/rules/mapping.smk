# Section 1 of the POC guideline: QC + combined-panel mapping.
# Ready to run once config/samples.tsv (read1/read2 paths) and
# data/references/panel_manifest.tsv (reference fasta paths) are filled in.
#
# Short-read (bwa-mem2, map_to_panel) and long-read (minimap2,
# map_longread_to_panel) isolates both land on the same output pattern
# (results/mapping/{sample}.sorted.bam) so every downstream rule (mosdepth,
# pav_call, ...) is aligner-agnostic. wildcard_constraints keep the two
# mapping rules from being ambiguous for the same {sample} value.

import re

_short_read_sample_pattern = "|".join(re.escape(s) for s in sample_ids) or "(?!)"
_long_read_sample_pattern = "|".join(re.escape(s) for s in sample_ids_longread) or "(?!)"

rule combine_reference_panel:
    input:
        fastas=references["fasta"].tolist(),
    output:
        "data/references/panel_combined.fa",
    shell:
        "cat {input.fastas} > {output}"


rule index_panel:
    input:
        "data/references/panel_combined.fa",
    output:
        multiext(
            "data/references/panel_combined.fa",
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
    log:
        "logs/mapping/index_panel.log",
    conda:
        "../../envs/mapping.yaml"
    shell:
        "bwa-mem2 index {input} > {log} 2>&1"


def sample_reads(wildcards):
    row = samples.set_index("sample_id").loc[wildcards.sample]
    return {"r1": row["read1"], "r2": row["read2"]}


rule fastp_qc:
    input:
        unpack(sample_reads),
    output:
        r1="results/mapping/{sample}.qc_R1.fastq.gz",
        r2="results/mapping/{sample}.qc_R2.fastq.gz",
        json="results/mapping/{sample}.fastp.json",
        html="results/mapping/{sample}.fastp.html",
    threads: 2
    log:
        "logs/mapping/{sample}_fastp.log",
    conda:
        "../../envs/mapping.yaml"
    shell:
        "fastp -i {input.r1} -I {input.r2} "
        "-o {output.r1} -O {output.r2} "
        "--json {output.json} --html {output.html} "
        "-w {threads} "
        "> {log} 2>&1"


rule map_to_panel:
    input:
        r1="results/mapping/{sample}.qc_R1.fastq.gz",
        r2="results/mapping/{sample}.qc_R2.fastq.gz",
        panel="data/references/panel_combined.fa",
        index=multiext(
            "data/references/panel_combined.fa",
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
    output:
        bam="results/mapping/{sample}.sorted.bam",
        bai="results/mapping/{sample}.sorted.bam.bai",
    threads: config["mapping"]["threads"]
    params:
        threads=config["mapping"]["threads"],
    log:
        "logs/mapping/{sample}_bwa.log",
    conda:
        "../../envs/mapping.yaml"
    wildcard_constraints:
        sample=_short_read_sample_pattern,
    shell:
        "bwa-mem2 mem -t {params.threads} {input.panel} {input.r1} {input.r2} 2> {log} "
        "| samtools sort -@4 -o {output.bam} - "
        "&& samtools index {output.bam}"


def longread_reads(wildcards):
    return samples_longread.set_index("sample_id").loc[wildcards.sample, "reads"]


def longread_minimap2_preset(wildcards):
    platform = samples_longread.set_index("sample_id").loc[wildcards.sample, "platform"]
    return {"nanopore": "map-ont", "pacbio": "map-pb"}.get(platform, "map-ont")


rule map_longread_to_panel:
    input:
        reads=longread_reads,
        panel="data/references/panel_combined.fa",
    output:
        bam="results/mapping/{sample}.sorted.bam",
        bai="results/mapping/{sample}.sorted.bam.bai",
    threads: config["mapping"]["threads"]
    params:
        threads=config["mapping"]["threads"],
        preset=longread_minimap2_preset,
    log:
        "logs/mapping/{sample}_minimap2.log",
    conda:
        "../../envs/mapping.yaml"
    wildcard_constraints:
        sample=_long_read_sample_pattern,
    shell:
        "minimap2 -t {params.threads} -ax {params.preset} {input.panel} {input.reads} 2> {log} "
        "| samtools sort -@4 -o {output.bam} - "
        "&& samtools index {output.bam}"
