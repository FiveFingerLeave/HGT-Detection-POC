# Assembles data/references/candidate_regions.bed (Section 1/3 of the POC
# guideline) from two sources:
# - "candidate": Starfish YR/Starship loci per reference (starfish.smk)
# - "core": core-marker windows per reference (generate_core_marker_regions.py)

rule core_marker_regions:
    input:
        fasta=lambda wc: references.set_index("reference_id").loc[wc.genome_id, "fasta"],
        exclude="results/starfish/starship_candidates.bed",
    output:
        "results/references/{genome_id}.core_regions.bed",
    log:
        "logs/references/{genome_id}_core_regions.log",
    conda:
        "../../envs/python.yaml"
    shell:
        "python workflow/scripts/generate_core_marker_regions.py "
        "--fasta {input.fasta} --genome-id {wildcards.genome_id} "
        "--exclude-bed {input.exclude} "
        "--output {output} "
        "> {log} 2>&1"


rule build_candidate_regions_bed:
    input:
        candidate="results/starfish/starship_candidates.bed",
        core=expand("results/references/{genome_id}.core_regions.bed", genome_id=reference_ids),
    output:
        "data/references/candidate_regions.bed",
    shell:
        "cat {input.candidate} {input.core} > {output}"
