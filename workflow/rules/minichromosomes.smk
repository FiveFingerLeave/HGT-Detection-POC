rule mask_repeats:
    input:
        fasta="results/normalized/{sample}.fna",
    output:
        masked="results/minichromosomes/{sample}.masked.fna",
    params:
        counts="results/minichromosomes/{sample}.windowmasker_counts.txt",
        conda_env=config["repeat_masking"]["conda_env"],
    log:
        "logs/minichromosomes/{sample}_windowmasker.log",
    shell:
        "conda run -n {params.conda_env} windowmasker -mk_counts "
        "-in {input.fasta} -out {params.counts} > {log} 2>&1 && "
        "conda run -n {params.conda_env} windowmasker -ustat {params.counts} "
        "-in {input.fasta} -out {output.masked} -outfmt fasta >> {log} 2>&1"


def contig_metrics_input(wildcards):
    inputs = {
        "fasta": f"results/normalized/{wildcards.sample}.fna",
        "id_map": f"results/normalized/{wildcards.sample}.id_map.tsv",
        "masked_fasta": f"results/minichromosomes/{wildcards.sample}.masked.fna",
    }
    if has_gff(wildcards.sample):
        inputs["gff"] = f"results/normalized/{wildcards.sample}.gff"
    if has_sequence_report(wildcards.sample):
        inputs["sequence_report"] = sample_sequence_report(wildcards)
    return inputs


def contig_metrics_args(wildcards):
    args = f" --masked-fasta results/minichromosomes/{wildcards.sample}.masked.fna"
    if has_gff(wildcards.sample):
        args += f" --gff results/normalized/{wildcards.sample}.gff"
    if has_sequence_report(wildcards.sample):
        args += (
            f" --sequence-report {sample_sequence_report(wildcards)}"
            f" --id-map results/normalized/{wildcards.sample}.id_map.tsv"
        )
    return args


rule calculate_contig_metrics:
    input:
        unpack(contig_metrics_input),
    output:
        "results/minichromosomes/{sample}.contig_metrics.tsv",
    params:
        extra_args=contig_metrics_args,
    log:
        "logs/minichromosomes/{sample}_contig_metrics.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/calculate_contig_metrics.py "
        "--fasta {input.fasta} --isolate-id {wildcards.sample} "
        "{params.extra_args} "
        "--output {output} "
        "> {log} 2>&1"


rule combine_contig_metrics:
    input:
        expand(
            "results/minichromosomes/{sample}.contig_metrics.tsv",
            sample=config["pilot_samples"],
        ),
    output:
        "results/minichromosomes/contig_metrics.unclassified.tsv",
    shell:
        "awk 'FNR==1 && NR!=1{{next}}{{print}}' {input} > {output}"


rule define_mchr_candidates:
    input:
        "results/minichromosomes/contig_metrics.unclassified.tsv",
    output:
        "results/minichromosomes/contig_metrics.tsv",
    params:
        core_min_length_bp=config["mchr"]["core_min_length_bp"],
        small_max_length_bp=config["mchr"]["small_max_length_bp"],
        gc_deviation_threshold=config["mchr"]["gc_deviation_threshold"],
        repeat_enrichment_threshold=config["mchr"]["repeat_enrichment_threshold"],
        gene_density_ratio_threshold=config["mchr"]["gene_density_ratio_threshold"],
    log:
        "logs/minichromosomes/define_mchr_candidates.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/define_mchr_candidates.py "
        "--input {input} --output {output} "
        "--core-min-length-bp {params.core_min_length_bp} "
        "--small-max-length-bp {params.small_max_length_bp} "
        "--gc-deviation-threshold {params.gc_deviation_threshold} "
        "--repeat-enrichment-threshold {params.repeat_enrichment_threshold} "
        "--gene-density-ratio-threshold {params.gene_density_ratio_threshold} "
        "> {log} 2>&1"
