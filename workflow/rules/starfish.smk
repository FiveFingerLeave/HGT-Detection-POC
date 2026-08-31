rule starfish_annotate_yr:
    input:
        fasta="results/normalized/{sample}.fna",
        ok="results/normalized/{sample}.validated.ok",
    output:
        gff="results/starfish/{sample}/{sample}_YR.filt.gff",
        fasta="results/starfish/{sample}/{sample}_YR.filt.fas",
    params:
        assembly_table="results/starfish/{sample}/{sample}_assembly.tsv",
        outdir="results/starfish/{sample}",
        tempdir="results/starfish/{sample}/tmp",
        profile=config["starfish"]["profile"],
        proteins=config["starfish"]["proteins"],
        idtag=config["starfish"]["idtag"],
        conda_env=config["starfish"]["conda_env"],
        separator=config["separator"],
    log:
        "logs/starfish/{sample}_annotate.log",
    threads: config["starfish"]["threads"]
    shell:
        "mkdir -p {params.outdir} {params.tempdir} && "
        "printf '%s\\t%s\\n' {wildcards.sample} {input.fasta} > {params.assembly_table} && "
        "conda run -n {params.conda_env} starfish annotate "
        "--assembly {params.assembly_table} "
        "--profile {params.profile} "
        "--proteins {params.proteins} "
        "--prefix {wildcards.sample}_YR "
        "--idtag {params.idtag} "
        "--outdir {params.outdir} "
        "--tempdir {params.tempdir} "
        "--threads {threads} "
        "--separator '{params.separator}' "
        "> {log} 2>&1"
