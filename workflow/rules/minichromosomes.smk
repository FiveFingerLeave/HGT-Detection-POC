def contig_metrics_input(wildcards):
    inputs = {
        "fasta": f"results/normalized/{wildcards.sample}.fna",
        "id_map": f"results/normalized/{wildcards.sample}.id_map.tsv",
    }
    if has_gff(wildcards.sample):
        inputs["gff"] = f"results/normalized/{wildcards.sample}.gff"
    if has_sequence_report(wildcards.sample):
        inputs["sequence_report"] = sample_sequence_report(wildcards)
    return inputs


def contig_metrics_args(wildcards):
    args = ""
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
        "results/minichromosomes/contig_metrics.tsv",
    shell:
        "awk 'FNR==1 && NR!=1{{next}}{{print}}' {input} > {output}"
