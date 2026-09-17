# Phase I (Section 6.2): uniform gene annotation across all 5 panel
# references (BRAKER3/Liftoff + InterProScan/eggNOG-mapper + OrthoFinder).
#
# Status: PARTIALLY implemented. BRAKER3 (ab-initio) remains blocked
# (needs a separately obtained GeneMark-ES/ET license file, see
# https://github.com/Gaius-Augustus/BRAKER - not a pure conda install).
# InterProScan/eggNOG-mapper (functional annotation) and OrthoFinder
# (Section 7.1) are also not yet implemented.
#
# Instead: Liftoff transfers the single available real annotation
# (GCA004346965_1, from NCBI) onto all 5 panel genomes (including itself,
# as a consistency check/renaming step). The document itself explicitly
# warns that lift-over alone can underestimate accessory genes, especially
# in non-syntenic/subtelomeric/repeat-rich regions - exactly where the four
# mini-/accessory-chromosome candidates found in repeats.smk are located.
# De-novo annotation of these regions remains an open point (see
# docs/decisions.md).


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
        # Dokument schreibt "-M msa -T iqtree" vor. Getestet und wieder
        # verworfen (siehe docs/decisions.md): mit nur 5 Spezies bleibt
        # die Orthogruppen-Zuordnung selbst (diamond+MCL) schnell (Minuten),
        # aber die anschliessende Gen-Baum-Inferenz laeuft iqtree3 MIT
        # ModelFinder PRO Orthogruppe einzeln - bei >11.000 Orthogruppen
        # und >15 Min je Baum waere das ein mehrtaegiger bis
        # mehrwoechiger Lauf. Fuer den in Abschnitt 7.1 tatsaechlich
        # benoetigten Output (die Orthogroups.tsv-Praesenz/Abwesenheits-
        # Matrix fuer die genbasierte PAV-Klassifikation) ist das nicht
        # noetig - "-M dendroblast" liefert dieselbe Orthogruppen-Matrix
        # (unveraendert durch diamond+MCL bestimmt) ohne die teure
        # MSA/Gen-Baum-Verfeinerungsschicht.
        "rm -rf {params.proteome_dir}/OrthoFinder {params.fixed_dir} && "
        "mkdir -p results/orthofinder && "
        "orthofinder -f {params.proteome_dir} -S diamond -M dendroblast "
        "-t {threads} -a {threads} > {log} 2>&1 && "
        "mv {params.proteome_dir}/OrthoFinder/Results_*/ {params.fixed_dir}"
