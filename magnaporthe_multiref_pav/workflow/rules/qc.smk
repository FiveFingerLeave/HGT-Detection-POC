# Phase I (Documentation/multireference_panel_pav_workflow.md, Section 6.1):
# assembly QC and header standardization of the 14 reference genomes.
# Per Section 17, immediately runnable (low effort).


rule rename_reference_headers:
    input:
        fasta=lambda wc: "data/references_raw/" + wc.genome_id + ".fna",
    output:
        fasta="data/references/{genome_id}.fa",
    params:
        name_map="results/panel/contig_name_map.tsv",
    conda:
        "../../envs/reporting.yaml"
    shell:
        "mkdir -p results/panel && "
        "python workflow/scripts/rename_fasta_headers.py "
        "--fasta {input.fasta} --genome-id {wildcards.genome_id} "
        "--output-fasta {output.fasta} --output-map {params.name_map}"


rule assembly_stats:
    input:
        expand("data/references/{genome_id}.fa", genome_id=genome_ids),
    output:
        "results/qc/assembly_stats.tsv",
    conda:
        "../../envs/core.yaml"
    shell:
        "seqkit stats -a -T {input} > {output}"


rule busco_reference:
    # A single BUSCO genome-mode run (metaeuk gene prediction against a
    # ~44 Mb genome) uses ~7 GB RSS - confirmed via an actual OOM-kill
    # when 2 ran in parallel on this 10 GB WSL2 VM (see docs/decisions.md).
    # threads: reserves the full core budget per job so Snakemake
    # serializes rather than oversubscribing, matching --cores 1.
    threads: config["threads_default"]
    input:
        "data/references/{genome_id}.fa",
    output:
        summary="results/qc/{genome_id}_busco_summary.txt",
    params:
        lineage=config["busco_lineage"],
        out_name="{genome_id}_busco",
        out_dir="results/qc",
        threads=config["threads_default"],
    log:
        "logs/qc/{genome_id}_busco.log",
    conda:
        "../../envs/core.yaml"
    shell:
        "cd {params.out_dir} && "
        "busco -f -i ../../{input} -l {params.lineage} -m genome "
        "-o {params.out_name} -c {params.threads} "
        "> ../../{log} 2>&1 && "
        "cd ../.. && "
        "find {params.out_dir}/{params.out_name} -name 'short_summary*.txt' "
        "-exec cp {{}} {output.summary} \\;"
