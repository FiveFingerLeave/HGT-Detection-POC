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


rule align_genome_self:
    input:
        fasta="results/normalized/{sample}.fna",
    output:
        paf="results/minichromosomes/{sample}.self_alignment.paf",
    params:
        conda_env=config["alignment"]["conda_env"],
        preset=config["alignment"]["minimap2_preset"],
        num_secondary=config["alignment"]["minimap2_num_secondary"],
    log:
        "logs/minichromosomes/{sample}_minimap2_self.log",
    threads: config["alignment"]["threads"]
    shell:
        "conda run -n {params.conda_env} minimap2 -t {threads} "
        "-x {params.preset} --secondary=yes -N {params.num_secondary} "
        "{input.fasta} {input.fasta} > {output.paf} 2> {log}"


rule compute_core_synteny:
    input:
        fasta="results/normalized/{sample}.fna",
        paf="results/minichromosomes/{sample}.self_alignment.paf",
    output:
        "results/minichromosomes/{sample}.core_synteny.tsv",
    params:
        core_min_length_bp=config["mchr"]["core_min_length_bp"],
    log:
        "logs/minichromosomes/{sample}_core_synteny.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/compute_core_synteny.py "
        "--paf {input.paf} --fasta {input.fasta} "
        "--core-min-length-bp {params.core_min_length_bp} "
        "--output {output} "
        "> {log} 2>&1"


def contig_metrics_input(wildcards):
    inputs = {
        "fasta": f"results/normalized/{wildcards.sample}.fna",
        "id_map": f"results/normalized/{wildcards.sample}.id_map.tsv",
        "masked_fasta": f"results/minichromosomes/{wildcards.sample}.masked.fna",
        "core_synteny": f"results/minichromosomes/{wildcards.sample}.core_synteny.tsv",
        "reference_synteny": f"results/minichromosomes/{wildcards.sample}.reference_synteny.tsv",
    }
    if has_gff(wildcards.sample):
        inputs["gff"] = f"results/normalized/{wildcards.sample}.gff"
    if has_sequence_report(wildcards.sample):
        inputs["sequence_report"] = sample_sequence_report(wildcards)
    return inputs


def contig_metrics_args(wildcards):
    args = (
        f" --masked-fasta results/minichromosomes/{wildcards.sample}.masked.fna"
        f" --telomere-window-bp {config['mchr']['telomere_window_bp']}"
        f" --telomere-min-repeats {config['mchr']['telomere_min_repeats']}"
        f" --core-synteny results/minichromosomes/{wildcards.sample}.core_synteny.tsv"
        f" --reference-synteny results/minichromosomes/{wildcards.sample}.reference_synteny.tsv"
    )
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
        core_synteny_veto_coverage=config["mchr"]["core_synteny_veto_coverage"],
        core_synteny_low_threshold=config["mchr"]["core_synteny_low_threshold"],
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
        "--core-synteny-veto-coverage {params.core_synteny_veto_coverage} "
        "--core-synteny-low-threshold {params.core_synteny_low_threshold} "
        "> {log} 2>&1"
