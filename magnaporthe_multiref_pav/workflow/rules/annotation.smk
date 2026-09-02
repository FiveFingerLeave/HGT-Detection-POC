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


# Section 7.1: Orthogruppen. gffread extrahiert die Proteinsequenzen aus
# den Liftoff-GFF3 (Voraussetzung fuer orthofinder -f), dann orthofinder
# selbst genau nach Dokument-Befehl (Threadzahlen an die 14 verfuegbaren
# Kerne angepasst statt der Dokument-Beispielwerte 32/16).


rule extract_proteome:
    # gffread -y markiert Stopcodons/nicht sauber uebersetzbare Codons
    # (z.B. an Exon-Grenzen mit Frame-Rest) mit "." - diamond akzeptiert
    # das nicht im Sequenzalphabet (siehe docs/decisions.md: diamond
    # v2.2.6 haengt sich bei so einem Zeichen sogar komplett auf statt
    # einen Fehler zu werfen). Ersetze "." (und vorsorglich "-") in
    # Sequenzzeilen (nicht Header) durch "X" (unbekannte Aminosaeure).
    input:
        fasta="data/references/{genome_id}.fa",
        gff="data/annotations/{genome_id}.gff3",
    output:
        "data/annotations/proteomes/{genome_id}.faa",
    log:
        "logs/annotation/{genome_id}_gffread.log",
    conda:
        "../../envs/orthofinder.yaml"
    shell:
        "mkdir -p data/annotations/proteomes && "
        "gffread -g {input.fasta} -y /dev/stdout {input.gff} 2> {log} "
        "| sed '/^>/!{{s/\\./X/g; s/-/X/g}}' > {output}"


rule orthofinder:
    # OrthoFinder's -o requires a not-yet-existing directory, which
    # conflicts with Snakemake auto-creating output parent dirs - so we
    # let it use its default OrthoFinder/Results_<timestamp>/ location
    # inside the proteome dir and move it to a fixed path afterwards.
    threads: config["threads_default"]
    input:
        expand("data/annotations/proteomes/{genome_id}.faa", genome_id=genome_ids),
    output:
        "results/orthofinder/Results/Orthogroups/Orthogroups.tsv",
    params:
        proteome_dir="data/annotations/proteomes",
        fixed_dir="results/orthofinder/Results",
    log:
        "logs/annotation/orthofinder.log",
    conda:
        "../../envs/orthofinder.yaml"
    shell:
        # Dokument schreibt "-T iqtree" vor; die installierte OrthoFinder-
        # Version (bundlet iqtree3 statt iqtree2) erwartet stattdessen den
        # Methodennamen "iqtree3" - inhaltlich dieselbe Wahl (IQ-TREE statt
        # FastTree/RAxML fuer Gen-Baeume), nur der CLI-Bezeichner hat sich
        # geaendert.
        "rm -rf {params.proteome_dir}/OrthoFinder {params.fixed_dir} && "
        "mkdir -p results/orthofinder && "
        "orthofinder -f {params.proteome_dir} -S diamond -M msa -T iqtree3 "
        "-t {threads} -a {threads} > {log} 2>&1 && "
        "mv {params.proteome_dir}/OrthoFinder/Results_*/ {params.fixed_dir}"
