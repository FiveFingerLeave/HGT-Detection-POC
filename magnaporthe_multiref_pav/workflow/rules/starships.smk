# Phase II (Section 7.4): Starship-like-Kandidaten (DUF3435-/Captain-
# Evidenz + Mindestgroesse + Cargo-Gene + fehlende Core-Syntenie).
#
# Schritt 1 (dieser Datei): `starfish annotate` sucht de-novo nach
# HMM-validierten Tyrosin-Rekombinase-("Captain")-Genen ueber alle 5
# Panel-Genome hinweg, in einem gemeinsamen Lauf (Multi-Genome-Assembly-
# TSV) - reproduziert den erfolgreichen Einzel-Isolat-Test aus
# Dokumentation/Starfish_bisherige_Schritte.md, jetzt fuer das ganze
# Panel. `-s '__'` passt den Separator an unsere bestehende
# {genome_id}__{contig}-Kopfzeilenkonvention an (Default waere '_',
# was bei genome_ids mit eigenem Unterstrich wie "GCA036493215_1"
# falsch parsen wuerde). `--gff` bindet die Liftoff-Genmodelle
# (annotation.smk) ein, damit spaetere Cargo-Gen-Analysen nicht nur auf
# den neu vorhergesagten YR-Genen basieren.
#
# Schritt 2 (noch offen): Verknuepfung der YR-Kandidaten mit
# Mindestgroesse, Cargo-Genen, Repeat-Kontext (repeats.smk) und
# fehlender Core-Syntenie (wga.smk) zur konservativen
# "starship_like"-Klassifikation gemaess Abschnitt 7.4.


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
    # Kombiniert die YR/Captain-Treffer mit Mindestgroesse (>=20kb,
    # thresholds.yaml: starship.min_region_length_bp), Repeat-Kontext
    # (repeats.smk) und Cargo-Genen (Liftoff-GFF3, annotation.smk) zur
    # konservativen Klassifikation nach Abschnitt 7.4. SyRI-Syntenie-
    # Kreuzreferenz (nur fuer die 3 kompatiblen Genome aus wga.smk)
    # ist noch nicht eingebaut - siehe docs/decisions.md.
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
    # Section 7.3: klassifiziert jedes 10kb-Fenster jedes Panel-Genoms
    # in einen Regionstyp, kombiniert aus Orthogruppen-Praevalenz (7.1),
    # Repeat-Dichte (6.3), Contig-Groesse/Gendichte (Mini-/Accessory-
    # Chromosom-Signal) und Starship-like-Kandidaten (7.4). SyRI-Syntenie
    # (nur 3/5 Genome, siehe wga.smk) ist noch NICHT eingebaut - siehe
    # docs/decisions.md fuer die Begruendung und als offener Folgeschritt.
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
