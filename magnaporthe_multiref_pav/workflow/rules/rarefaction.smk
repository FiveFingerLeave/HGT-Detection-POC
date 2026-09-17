# Phase IX (Section 14): Reference rarefaction (14.1/14.2) and
# test-isolate rarefaction/novelty check (14.3, reduced scope - see
# below and docs/decisions.md).
#
# 14.1/14.2: The document specifies random k-out-of-14 combinations
# (1000 permutations), because C(14,k) would be too large for exhaustive
# enumeration. Our panel, deliberately reduced to 5 host representatives,
# instead allows EXHAUSTIVE enumeration of all C(5,k) combinations
# (at most 10) - stricter than the document's own sampling procedure,
# not a weakening.
#
# 14.3: Full implementation (unmapped reads -> local assembly -> panel
# lookup -> new candidate regions) requires a long-read assembler
# (e.g. Flye), which is not yet installed. Only the first, cheap
# part (unmapped-read extraction + basic statistics) is implemented
# here, as an approximation for the question "how much isolate sequence
# does the panel not explain" - the actual contig assembly/novelty
# classification is NOT implemented.


rule rarefaction_reference_panel:
    input:
        manifest="results/panel/panel_contig_manifest.tsv",
    output:
        table="results/rarefaction/rarefaction_reference_panel.tsv",
        saturation="results/rarefaction/saturation_summary.tsv",
    params:
        genome_ids=genome_ids,
    log:
        "logs/rarefaction/reference_panel.log",
    conda:
        "../../envs/core.yaml"
    shell:
        "mkdir -p results/rarefaction logs/rarefaction && "
        "python3 workflow/scripts/rarefaction_reference_panel.py "
        "--manifest {input.manifest} "
        "--genome-ids {params.genome_ids} "
        "--out {output.table} "
        "--out-saturation {output.saturation} "
        "> {log} 2>&1"


rule extract_unmapped_reads:
    # Section 14.3, step 1 (partial implementation): reads that could
    # not be assigned to the panel at all.
    threads: config["threads_default"]
    input:
        bam="results/mapping/{sample_id}.panel.bam",
        bai="results/mapping/{sample_id}.panel.bam.bai",
    output:
        fastq="results/rarefaction/{sample_id}_unmapped.fastq.gz",
    conda:
        "../../envs/core.yaml"
    shell:
        "mkdir -p results/rarefaction && "
        "samtools view -@ {threads} -b -f 4 {input.bam} "
        "| samtools fastq -@ {threads} - 2>/dev/null | gzip > {output.fastq}"


rule unmapped_read_stats:
    input:
        fastq="results/rarefaction/{sample_id}_unmapped.fastq.gz",
    output:
        stats="results/rarefaction/{sample_id}_unmapped_stats.tsv",
    conda:
        "../../envs/core.yaml"
    shell:
        "seqkit stats -a {input.fastq} > {output.stats}"


rule rarefaction_all:
    input:
        "results/rarefaction/rarefaction_reference_panel.tsv",
        "results/rarefaction/saturation_summary.tsv",
        expand("results/rarefaction/{sample_id}_unmapped_stats.tsv", sample_id=mappable_sample_ids),
