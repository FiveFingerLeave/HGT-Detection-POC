# Section 2 of the POC guideline: structure-based Starship detection.
# Runs over every genome (reference panel + any isolate assemblies) listed
# in data/references/panel_manifest.tsv - add isolate assemblies there too
# once available, or extend with a second manifest for data/isolates_poc/.

rule starfish_annotate:
    input:
        assembly=lambda wc: references.set_index("reference_id").loc[wc.genome_id, "fasta"],
        gff=lambda wc: references.set_index("reference_id").loc[wc.genome_id, "gff"],
    output:
        directory("results/starfish/{genome_id}"),
    params:
        conda_env=config["starfish"]["conda_env"],
        idtag=config["starfish"]["idtag"],
        threads=config["starfish"]["threads"],
    log:
        "logs/starfish/{genome_id}_annotate.log",
    shell:
        "conda run -n {params.conda_env} starfish annotate -T {params.threads} "
        "-x {wildcards.genome_id} -i {params.idtag} "
        "-a {input.assembly} -g {input.gff} "
        "-o {output} "
        "> {log} 2>&1"


rule starfish_insert:
    input:
        "results/starfish/{genome_id}",
    output:
        directory("results/starfish/{genome_id}.inserts"),
    params:
        conda_env=config["starfish"]["conda_env"],
    log:
        "logs/starfish/{genome_id}_insert.log",
    shell:
        "conda run -n {params.conda_env} starfish insert -x {wildcards.genome_id} "
        "-a {input}/{wildcards.genome_id}.starships.bed "
        "-o {output} "
        "> {log} 2>&1"


rule stargraph_build:
    input:
        assemblies=references["fasta"].tolist(),
        starships=expand(
            "results/starfish/{genome_id}", genome_id=reference_ids
        ),
    output:
        directory("results/starfish/pangenome_graph"),
    params:
        conda_env=config["starfish"]["conda_env"],
    log:
        "logs/starfish/stargraph_build.log",
    shell:
        # Stargraph is not on conda/bioconda - install from
        # https://github.com/egluckthaler/stargraph into {params.conda_env} first.
        "conda run -n {params.conda_env} stargraph build "
        "--assemblies {input.assemblies} "
        "--starships {input.starships}/*.bed "
        "-o {output} "
        "> {log} 2>&1"
