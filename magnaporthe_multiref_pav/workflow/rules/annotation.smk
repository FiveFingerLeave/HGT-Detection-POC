# Phase I (Section 6.2): einheitliche Genannotation ueber alle 5
# Panel-Referenzen (BRAKER3/Liftoff + InterProScan/eggNOG-mapper +
# OrthoFinder).
#
# Status: TEILWEISE implementiert. BRAKER3 (ab-initio) bleibt blockiert
# (braucht eine separat zu beschaffende GeneMark-ES/ET-Lizenzdatei, siehe
# https://github.com/Gaius-Augustus/BRAKER - kein reiner conda-Install).
# InterProScan/eggNOG-mapper (Funktionsannotation) und OrthoFinder
# (Abschnitt 7.1) ebenfalls noch nicht implementiert.
#
# Stattdessen: Liftoff ueberträgt die einzige vorhandene echte Annotation
# (GCA004346965_1, aus NCBI) auf alle 5 Panel-Genome (inkl. sich selbst,
# als Konsistenzpruefung/Umbenennungsschritt). Das Dokument selbst warnt
# ausdruecklich, dass Lift-over allein akzessorische Gene unterschaetzen
# kann, besonders in nichtsyntenischen/subtelomerischen/repeat-reichen
# Bereichen - genau dort, wo die vier in repeats.smk gefundenen
# Mini-/Accessory-Chromosom-Kandidaten liegen. De-novo-Annotation dieser
# Bereiche bleibt ein offener Punkt (siehe docs/decisions.md).


rule liftoff_annotation:
    threads: config["threads_default"]
    input:
        target_fasta="data/references/{genome_id}.fa",
        ref_fasta_raw="data/references_raw/GCA004346965_1.fna",
        ref_gff="data/annotations_raw/GCA004346965_1.gff3",
    output:
        gff="data/annotations/{genome_id}.gff3",
        unmapped="results/annotation/{genome_id}_unmapped_features.txt",
    params:
        outdir="results/annotation",
        intermediate="results/annotation/{genome_id}_intermediate",
    log:
        "logs/annotation/{genome_id}_liftoff.log",
    conda:
        "../../envs/liftoff.yaml"
    shell:
        "mkdir -p {params.outdir} data/annotations && "
        "liftoff -g {input.ref_gff} -o {output.gff} -u {output.unmapped} "
        "-dir {params.intermediate} "
        "-p {threads} {input.target_fasta} {input.ref_fasta_raw} "
        "> {log} 2>&1"


rule annotation_all:
    input:
        expand("data/annotations/{genome_id}.gff3", genome_id=genome_ids),
