# Section 2 of the POC guideline: structure-based Starship/YR detection
# (starfish annotate) on every reference-panel genome. The resulting
# filtered YR/Starship candidate loci become the region set tested for
# coverage-based presence/absence in Section 1 (cross-validation between
# methods a+b, as the guideline specifies).
#
# NOTE: reference FASTA headers are NOT genome-ID-prefixed (no
# normalization step in this rebuild). starfish annotate will warn
# ("... is being parsed into <2 components using separator ...") but still
# completes and finds real YR candidates - confirmed empirically on this
# same pilot isolate earlier in the project (see docs/decisions.md). This
# is fine for the current scope (annotate -> candidate coordinates ->
# candidate_regions.bed): those coordinates use the same un-prefixed
# contig names as data/references/panel_combined.fa, so they still line
# up correctly for PAV calling.

rule write_starfish_assembly_table:
    input:
        fasta=lambda wc: references.set_index("reference_id").loc[wc.genome_id, "fasta"],
    output:
        "results/starfish/{genome_id}/{genome_id}_assembly.tsv",
    shell:
        "mkdir -p $(dirname {output}) && printf '%s\\t%s\\n' {wildcards.genome_id} {input.fasta} > {output}"


def starfish_gff(genome_id: str) -> str:
    gff = references.set_index("reference_id").loc[genome_id, "gff"]
    return gff if isinstance(gff, str) and gff.strip() else ""


def starfish_annotate_input(wildcards):
    inputs = {
        "assembly_table": f"results/starfish/{wildcards.genome_id}/{wildcards.genome_id}_assembly.tsv",
        "assembly": references.set_index("reference_id").loc[wildcards.genome_id, "fasta"],
    }
    gff = starfish_gff(wildcards.genome_id)
    if gff:
        inputs["gff"] = gff
    return inputs


def starfish_gff_args(wildcards):
    gff = starfish_gff(wildcards.genome_id)
    if not gff:
        return ""
    gff_table = f"results/starfish/{wildcards.genome_id}/{wildcards.genome_id}_gff.tsv"
    return f"-g {gff_table}"


rule starfish_annotate:
    input:
        unpack(starfish_annotate_input),
    output:
        gff="results/starfish/{genome_id}/{genome_id}_YR.filt.gff",
        fasta="results/starfish/{genome_id}/{genome_id}_YR.filt.fas",
    params:
        conda_env=config["starfish"]["conda_env"],
        profile=config["starfish"]["profile"],
        proteins=config["starfish"]["proteins"],
        idtag=config["starfish"]["idtag"],
        separator=config["starfish"]["separator"],
        threads=config["starfish"]["threads"],
        outdir="results/starfish/{genome_id}",
        tempdir="results/starfish/{genome_id}/tmp",
        gff_table="results/starfish/{genome_id}/{genome_id}_gff.tsv",
        gff_source=lambda wc: starfish_gff(wc.genome_id),
        gff_args=starfish_gff_args,
    log:
        "logs/starfish/{genome_id}_annotate.log",
    shell:
        "mkdir -p {params.tempdir} && "
        "if [ -n '{params.gff_source}' ]; then "
        "printf '%s\\t%s\\n' {wildcards.genome_id} {params.gff_source} > {params.gff_table}; "
        "fi && "
        "conda run -n {params.conda_env} starfish annotate "
        "--assembly {input.assembly_table} "
        "--profile {params.profile} "
        "--proteins {params.proteins} "
        "--prefix {wildcards.genome_id}_YR "
        "--idtag {params.idtag} "
        "--outdir {params.outdir} "
        "--tempdir {params.tempdir} "
        "--threads {params.threads} "
        "--separator '{params.separator}' "
        "{params.gff_args} "
        "> {log} 2>&1"


rule starship_candidates_to_bed:
    # Turn every reference's filtered YR/Starship candidate GFF into BED
    # regions for candidate_regions.bed (region_class=candidate).
    input:
        expand(
            "results/starfish/{genome_id}/{genome_id}_YR.filt.gff", genome_id=reference_ids
        ),
    output:
        "results/starfish/starship_candidates.bed",
    params:
        gff_args=" ".join(
            f"--gff {gid}=results/starfish/{gid}/{gid}_YR.filt.gff" for gid in reference_ids
        ),
    conda:
        "../../envs/python.yaml"
    shell:
        "python workflow/scripts/gff_to_bed.py "
        "{params.gff_args} "
        "--region-class candidate "
        "--id-prefix starship "
        "--output {output}"
