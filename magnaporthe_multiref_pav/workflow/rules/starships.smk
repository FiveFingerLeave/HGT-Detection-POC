# Phase II (Section 7.4): Starship-like candidates (DUF3435/Captain
# evidence + minimum size + cargo genes + missing core synteny).
#
# Step 1 (this file): `starfish annotate` searches de novo for
# HMM-validated tyrosine recombinase ("Captain") genes across all 5
# panel genomes, in one combined run (multi-genome assembly
# TSV) - reproduces the successful single-isolate test from
# Documentation/Starfish_Progress_Log.md, now for the whole
# panel. `-s '__'` adapts the separator to our existing
# {genome_id}__{contig} header convention (the default would be '_',
# which would parse incorrectly for genome_ids containing their own
# underscore, like "GCA036493215_1"). `--gff` incorporates the Liftoff
# gene models (annotation.smk) so that later cargo-gene analyses are
# not based solely on the newly predicted YR genes.
#
# Step 2 (still open): linking the YR candidates with minimum size,
# cargo genes, repeat context (repeats.smk) and missing core
# synteny (wga.smk) for the conservative "starship_like"
# classification per Section 7.4.


rule starfish_input_lists:
    output:
        assembly_tsv="results/starships/assembly_list.tsv",
        gff_tsv="results/starships/gff_list.tsv",
    run:
        import os

        os.makedirs("results/starships", exist_ok=True)
        with open(output.assembly_tsv, "w") as f:
            for g in genome_ids:
                f.write(f"{g}\tdata/references/{g}.fa\n")
        with open(output.gff_tsv, "w") as f:
            for g in genome_ids:
                f.write(f"{g}\tdata/annotations/{g}.gff3\n")


rule starfish_annotate_yr:
    threads: config["threads_default"]
    input:
        assembly_tsv="results/starships/assembly_list.tsv",
        gff_tsv="results/starships/gff_list.tsv",
        fastas=expand("data/references/{genome_id}.fa", genome_id=genome_ids),
        gffs=expand("data/annotations/{genome_id}.gff3", genome_id=genome_ids),
    output:
        gff="results/starships/panel_YR.filt.gff",
        fasta="results/starships/panel_YR.filt.fas",
    params:
        outdir="results/starships",
        tempdir="results/starships/tmp",
    log:
        "logs/starships/annotate_YR.log",
    conda:
        "../../envs/starships.yaml"
    shell:
        "mkdir -p {params.outdir} {params.tempdir} && "
        "starfish annotate "
        "--assembly {input.assembly_tsv} "
        "--gff {input.gff_tsv} "
        "--profile $CONDA_PREFIX/db/YRsuperfams.p1-512.hmm "
        "--proteins $CONDA_PREFIX/db/YRsuperfamRefs.faa "
        "--prefix panel_YR --idtag YR -s '__' "
        "--outdir {params.outdir} --tempdir {params.tempdir} "
        "--threads {threads} > {log} 2>&1"


rule classify_starship_candidates:
    # Combines the YR/Captain hits with minimum size (>=20kb,
    # thresholds.yaml: starship.min_region_length_bp), repeat context
    # (repeats.smk) and cargo genes (Liftoff GFF3, annotation.smk) for
    # the conservative classification per Section 7.4. SyRI synteny
    # cross-referencing (only for the 3 compatible genomes from wga.smk)
    # is not yet implemented - see docs/decisions.md.
    input:
        yr_gff="results/starships/panel_YR.filt.gff",
        fais=expand("data/references/{genome_id}.fa.fai", genome_id=genome_ids),
        repeat_windows=expand("results/repeats/{genome_id}_repeat_windows.bed", genome_id=genome_ids),
        repeat_per_contig=expand("results/repeats/{genome_id}_repeat_per_contig.tsv", genome_id=genome_ids),
        gffs=expand("data/annotations/{genome_id}.gff3", genome_id=genome_ids),
    output:
        "results/starships/starship_like_candidates.tsv",
    params:
        genome_ids=genome_ids,
        min_region_bp=thresholds["starship"]["min_region_length_bp"],
    log:
        "logs/starships/classify.log",
    conda:
        "../../envs/core.yaml"
    shell:
        "python3 workflow/scripts/classify_starship_candidates.py "
        "--yr-gff {input.yr_gff} "
        "--genome-ids {params.genome_ids} "
        "--fai-template 'data/references/{{genome_id}}.fa.fai' "
        "--repeat-windows-template 'results/repeats/{{genome_id}}_repeat_windows.bed' "
        "--repeat-per-contig-template 'results/repeats/{{genome_id}}_repeat_per_contig.tsv' "
        "--gff-template 'data/annotations/{{genome_id}}.gff3' "
        "--min-region-bp {params.min_region_bp} "
        "--out {output} > {log} 2>&1"


rule classify_panel_regions:
    # Section 7.3: classifies every 10kb window of every panel genome
    # into a region type, combining orthogroup prevalence (7.1),
    # repeat density (6.3), contig size/gene density (mini-/accessory-
    # chromosome signal) and starship-like candidates (7.4). SyRI synteny
    # (only 3/5 genomes, see wga.smk) is NOT yet implemented - see
    # docs/decisions.md for the rationale and as an open follow-up step.
    input:
        fais=expand("data/references/{genome_id}.fa.fai", genome_id=genome_ids),
        repeat_windows=expand("results/repeats/{genome_id}_repeat_windows.bed", genome_id=genome_ids),
        repeat_per_contig=expand("results/repeats/{genome_id}_repeat_per_contig.tsv", genome_id=genome_ids),
        gffs=expand("data/annotations/{genome_id}.gff3", genome_id=genome_ids),
        orthogroups="results/orthofinder/Results/Orthogroups/Orthogroups.tsv",
        starship_candidates="results/starships/starship_like_candidates.tsv",
    output:
        "results/panel/panel_regions.bed",
    params:
        genome_ids=genome_ids,
        strict_core=thresholds["panel"]["strict_core_fraction"],
        soft_core=thresholds["panel"]["soft_core_fraction"],
        shell_min=thresholds["panel"]["shell_min_fraction"],
    log:
        "logs/panel/classify_regions.log",
    conda:
        "../../envs/core.yaml"
    shell:
        "mkdir -p results/panel logs/panel && "
        "python3 workflow/scripts/classify_panel_regions.py "
        "--genome-ids {params.genome_ids} "
        "--fai-template 'data/references/{{genome_id}}.fa.fai' "
        "--repeat-windows-template 'results/repeats/{{genome_id}}_repeat_windows.bed' "
        "--repeat-per-contig-template 'results/repeats/{{genome_id}}_repeat_per_contig.tsv' "
        "--gff-template 'data/annotations/{{genome_id}}.gff3' "
        "--orthogroups {input.orthogroups} "
        "--starship-candidates {input.starship_candidates} "
        "--strict-core-fraction {params.strict_core} "
        "--soft-core-fraction {params.soft_core} "
        "--shell-min-fraction {params.shell_min} "
        "--out {output} > {log} 2>&1"
