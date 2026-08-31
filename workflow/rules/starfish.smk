def starfish_annotate_input(wildcards):
    inputs = {"fasta": f"results/normalized/{wildcards.sample}.fna"}
    if has_gff(wildcards.sample):
        inputs["validated"] = f"results/normalized/{wildcards.sample}.validated.ok"
        inputs["gff"] = f"results/normalized/{wildcards.sample}.gff"
    return inputs


def starfish_gff_source(wildcards):
    return f"results/normalized/{wildcards.sample}.gff" if has_gff(wildcards.sample) else ""


def starfish_gff_args(wildcards):
    if not has_gff(wildcards.sample):
        return ""
    return f"--gff results/starfish/{wildcards.sample}/{wildcards.sample}_gff.tsv"


rule starfish_annotate_yr:
    input:
        unpack(starfish_annotate_input),
    output:
        gff="results/starfish/{sample}/{sample}_YR.filt.gff",
        fasta="results/starfish/{sample}/{sample}_YR.filt.fas",
    params:
        assembly_table="results/starfish/{sample}/{sample}_assembly.tsv",
        gff_table="results/starfish/{sample}/{sample}_gff.tsv",
        gff_source=starfish_gff_source,
        gff_args=starfish_gff_args,
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
        "if [ -n '{params.gff_source}' ]; then "
        "printf '%s\\t%s\\n' {wildcards.sample} {params.gff_source} > {params.gff_table}; "
        "fi && "
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
        "{params.gff_args} "
        "> {log} 2>&1"
