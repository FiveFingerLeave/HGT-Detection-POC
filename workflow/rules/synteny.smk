# Synteny of Starship/accessory candidate regions across host lineages:
# pairwise whole-genome alignment (minimap2) between every pair of
# reference-panel genomes, then classify each candidate/core region by
# whether it is syntenic (aligns in the other genome), a lineage-specific
# insertion at an otherwise-conserved locus ("syntenic_empty_site" - the
# classic mobile-element insertion-polymorphism signature: flanks align,
# cargo does not), or uninformative (neither region nor flanks align).

import re

_reference_id_pattern = "|".join(re.escape(r) for r in reference_ids) or "(?!)"


def synteny_pairs():
    return [(g, o) for g in reference_ids for o in reference_ids if g != o]


rule pairwise_genome_alignment:
    input:
        query=lambda wc: references.set_index("reference_id").loc[wc.genome, "fasta"],
        target=lambda wc: references.set_index("reference_id").loc[wc.other_genome, "fasta"],
    output:
        "results/synteny/{genome}_vs_{other_genome}.paf",
    wildcard_constraints:
        genome=_reference_id_pattern,
        other_genome=_reference_id_pattern,
    params:
        threads=config["mapping"]["threads"],
    log:
        "logs/synteny/{genome}_vs_{other_genome}_minimap2.log",
    conda:
        "../../envs/mapping.yaml"
    shell:
        "minimap2 -x asm20 -t {params.threads} {input.target} {input.query} "
        "> {output} 2> {log}"


rule synteny_overlap:
    input:
        regions="data/references/candidate_regions.bed",
        fasta=lambda wc: references.set_index("reference_id").loc[wc.genome, "fasta"],
        paf="results/synteny/{genome}_vs_{other_genome}.paf",
    output:
        "results/synteny/{genome}_vs_{other_genome}.synteny.tsv",
    wildcard_constraints:
        genome=_reference_id_pattern,
        other_genome=_reference_id_pattern,
    log:
        "logs/synteny/{genome}_vs_{other_genome}_overlap.log",
    conda:
        "../../envs/python.yaml"
    shell:
        "python workflow/scripts/synteny_overlap.py "
        "--regions {input.regions} --fasta {input.fasta} --paf {input.paf} "
        "--genome {wildcards.genome} --other-genome {wildcards.other_genome} "
        "--output {output} "
        "> {log} 2>&1"


rule combine_synteny:
    input:
        [
            f"results/synteny/{g}_vs_{o}.synteny.tsv"
            for g, o in synteny_pairs()
        ],
    output:
        "results/synteny/candidate_synteny.tsv",
    run:
        import pandas as pd

        pd.concat(
            [pd.read_csv(f, sep="\t") for f in input], ignore_index=True
        ).to_csv(output[0], sep="\t", index=False)
