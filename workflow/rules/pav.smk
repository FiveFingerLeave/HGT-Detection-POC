# Section 1 of the POC guideline: coverage/breadth -> presence-absence
# calls per candidate region. Requires data/references/candidate_regions.bed
# (regions to test - e.g. known mChr/Starship loci from the reference panel).

rule mosdepth_coverage:
    input:
        bam="results/mapping/{sample}.sorted.bam",
        bai="results/mapping/{sample}.sorted.bam.bai",
        regions="data/references/candidate_regions.bed",
    output:
        regions="results/pav_calls/{sample}.regions.bed.gz",
        thresholds="results/pav_calls/{sample}.thresholds.bed.gz",
    params:
        prefix="results/pav_calls/{sample}",
        present_depth=config["pav"]["present_depth"],
    log:
        "logs/pav_calls/{sample}_mosdepth.log",
    conda:
        "../../envs/mapping.yaml"
    shell:
        "mosdepth --by {input.regions} -t 4 "
        "--thresholds {params.present_depth} "
        "{params.prefix} {input.bam} > {log} 2>&1"


rule pav_call:
    input:
        regions="data/references/candidate_regions.bed",
        mosdepth_regions="results/pav_calls/{sample}.regions.bed.gz",
        mosdepth_thresholds="results/pav_calls/{sample}.thresholds.bed.gz",
    output:
        "results/pav_calls/{sample}.pav.tsv",
    params:
        sample="{sample}",
        present_breadth=config["pav"]["present_breadth"],
        present_depth=config["pav"]["present_depth"],
        absent_breadth=config["pav"]["absent_breadth"],
    log:
        "logs/pav_calls/{sample}_pav_call.log",
    conda:
        "../../envs/python.yaml"
    shell:
        "python workflow/scripts/pav_call.py "
        "--sample-id {params.sample} "
        "--regions {input.regions} "
        "--mosdepth-regions {input.mosdepth_regions} "
        "--mosdepth-thresholds {input.mosdepth_thresholds} "
        "--present-breadth {params.present_breadth} "
        "--present-depth {params.present_depth} "
        "--absent-breadth {params.absent_breadth} "
        "--out {output} "
        "> {log} 2>&1"


rule combine_pav_calls:
    input:
        expand("results/pav_calls/{sample}.pav.tsv", sample=sample_ids),
    output:
        "results/pav_calls/candidate_table.tsv",
    shell:
        "awk 'FNR==1 && NR!=1{{next}}{{print}}' {input} > {output}"
