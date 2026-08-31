rule reference_core_contigs:
    output:
        "results/references/{reference_id}.core_contigs.txt",
    params:
        fasta=lambda wc: reference_by_id[wc.reference_id]["fasta"],
        sequence_report=lambda wc: reference_sequence_report(wc.reference_id),
        core_min_length_bp=config["mchr"]["core_min_length_bp"],
        exclude=lambda wc: ",".join(known_accessory_by_reference.get(wc.reference_id, [])),
    log:
        "logs/references/{reference_id}_core_contigs.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/list_reference_core_contigs.py "
        "--fasta {params.fasta} "
        "--sequence-report {params.sequence_report} "
        "--core-min-length-bp {params.core_min_length_bp} "
        "--exclude '{params.exclude}' "
        "--output {output} "
        "> {log} 2>&1"


rule align_to_reference:
    input:
        query="results/normalized/{sample}.fna",
    output:
        paf="results/minichromosomes/{sample}.vs_{reference_id}.paf",
    params:
        reference_fasta=lambda wc: reference_by_id[wc.reference_id]["fasta"],
        conda_env=config["alignment"]["conda_env"],
        preset=config["alignment"]["reference_minimap2_preset"],
    log:
        "logs/minichromosomes/{sample}_vs_{reference_id}_minimap2.log",
    threads: config["alignment"]["threads"]
    shell:
        "conda run -n {params.conda_env} minimap2 -t {threads} "
        "-x {params.preset} --secondary=no "
        "{params.reference_fasta} {input.query} > {output.paf} 2> {log}"


def compute_reference_synteny_input(wildcards):
    ref_ids = structural_reference_ids_for_sample(wildcards.sample)
    return {
        "fasta": f"results/normalized/{wildcards.sample}.fna",
        "pafs": [
            f"results/minichromosomes/{wildcards.sample}.vs_{rid}.paf" for rid in ref_ids
        ],
        "core_contig_lists": [
            f"results/references/{rid}.core_contigs.txt" for rid in ref_ids
        ],
    }


rule compute_reference_synteny:
    input:
        unpack(compute_reference_synteny_input),
    output:
        "results/minichromosomes/{sample}.reference_synteny.tsv",
    params:
        min_coverage=config["alignment"]["min_coverage_per_reference"],
        min_mapq=config["alignment"]["min_mapq"],
        consensus_min_hits=consensus_min_hits,
        paf_args=lambda wc, input: " ".join(f"--reference-paf {p}" for p in input.pafs),
        core_args=lambda wc, input: " ".join(
            f"--reference-core-contigs {c}" for c in input.core_contig_lists
        ),
    log:
        "logs/minichromosomes/{sample}_reference_synteny.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python workflow/scripts/compute_reference_synteny.py "
        "--query-fasta {input.fasta} "
        "{params.paf_args} "
        "{params.core_args} "
        "--min-coverage-per-reference {params.min_coverage} "
        "--min-mapq {params.min_mapq} "
        "--consensus-min-hits {params.consensus_min_hits} "
        "--output {output} "
        "> {log} 2>&1"
