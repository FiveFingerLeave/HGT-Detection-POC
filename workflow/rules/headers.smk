rule normalize_fasta_headers:
    input:
        fasta=sample_fasta,
    output:
        fasta="results/normalized/{sample}.fna",
        mapping="results/normalized/{sample}.id_map.tsv",
    log:
        "logs/headers/{sample}_fasta_headers.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/normalize_fasta_headers.py "
        "--fasta {input.fasta} --isolate-id {wildcards.sample} "
        "--output-fasta {output.fasta} --output-mapping {output.mapping} "
        "--separator '{config[separator]}' "
        "> {log} 2>&1"


rule normalize_gff_seqids:
    input:
        gff=sample_gff,
        mapping="results/normalized/{sample}.id_map.tsv",
    output:
        gff="results/normalized/{sample}.gff",
    log:
        "logs/headers/{sample}_gff_seqids.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/normalize_gff_seqids.py "
        "--gff {input.gff} --mapping {input.mapping} "
        "--output-gff {output.gff} "
        "> {log} 2>&1"


rule validate_fasta_gff_ids:
    input:
        fasta="results/normalized/{sample}.fna",
        gff="results/normalized/{sample}.gff",
    output:
        touch("results/normalized/{sample}.validated.ok"),
    log:
        "logs/headers/{sample}_validate.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/validate_fasta_gff_ids.py "
        "--fasta {input.fasta} --gff {input.gff} "
        "> {log} 2>&1"
