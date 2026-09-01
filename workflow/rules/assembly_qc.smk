rule run_quast:
    input:
        fastas=list(qc_genomes.values()),
    output:
        report="results/qc/quast/report.tsv",
    params:
        conda_env=config["qc"]["conda_env"],
        labels=",".join(qc_genomes.keys()),
        threads=config["qc"]["threads"],
        outdir="results/qc/quast",
    log:
        "logs/qc/quast.log",
    shell:
        "conda run -n {params.conda_env} quast.py "
        "{input.fastas} "
        "--labels {params.labels} "
        "--threads {params.threads} "
        "-o {params.outdir} "
        "> {log} 2>&1"


rule run_busco:
    input:
        fasta=lambda wc: qc_genomes[wc.genome_id],
    output:
        summary="results/qc/busco/{genome_id}/short_summary_{genome_id}.txt",
    params:
        conda_env=config["qc"]["conda_env"],
        lineage=config["qc"]["busco_lineage"],
        threads=config["qc"]["threads"],
        out_path="results/qc/busco",
    log:
        "logs/qc/busco_{genome_id}.log",
    shell:
        "conda run -n {params.conda_env} busco "
        "-i {input.fasta} "
        "-o {wildcards.genome_id} "
        "--out_path {params.out_path} "
        "-l {params.lineage} "
        "-m genome "
        "-c {params.threads} "
        "-f "
        "> {log} 2>&1 && "
        "cp {params.out_path}/{wildcards.genome_id}/short_summary.*.{wildcards.genome_id}.txt "
        "{output.summary}"


rule summarize_assembly_qc:
    input:
        quast_report="results/qc/quast/report.tsv",
        busco_summaries=expand(
            "results/qc/busco/{genome_id}/short_summary_{genome_id}.txt",
            genome_id=qc_genomes.keys(),
        ),
    output:
        "results/qc/assembly_qc_summary.tsv",
    params:
        busco_args=lambda wc, input: " ".join(
            f"--busco-summary {genome_id}={path}"
            for genome_id, path in zip(qc_genomes.keys(), input.busco_summaries)
        ),
        tier_args=" ".join(f"--tier {gid}={tier}" for gid, tier in qc_tiers.items()),
    log:
        "logs/qc/summarize_assembly_qc.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/summarize_assembly_qc.py "
        "--quast-report {input.quast_report} "
        "{params.busco_args} "
        "{params.tier_args} "
        "--output {output} "
        "> {log} 2>&1"
