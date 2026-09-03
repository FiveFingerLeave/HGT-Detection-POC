# Phase IX (Section 14): Referenz-Rarefaction (14.1/14.2) und
# Testisolat-Rarefaction/Novelty-Check (14.3, reduzierter Umfang - siehe
# unten und docs/decisions.md).
#
# 14.1/14.2: Das Dokument sieht zufaellige k-aus-14-Kombinationen vor
# (1000 Permutationen), weil C(14,k) zu gross fuer erschoepfende
# Aufzaehlung waere. Unser bewusst auf 5 Host-Repraesentanten reduziertes
# Panel erlaubt dagegen ERSCHOEPFENDE Aufzaehlung aller C(5,k)-Kombinationen
# (maximal 10) - strenger als das Dokument-eigene Sampling-Verfahren,
# keine Abschwaechung.
#
# 14.3: Volle Umsetzung (unmapped Reads -> lokale Assembly -> Panel-
# Ruecksuche -> neue Kandidatenregionen) braucht einen Long-Read-Assembler
# (z.B. Flye), der bisher nicht installiert ist. Umgesetzt ist hier nur
# der erste, guenstige Teil (unmapped-Read-Extraktion + Basisstatistik)
# als Naeherung fuer die Frage "wie viel Isolat-Sequenz erklaert das
# Panel nicht" - die eigentliche Contig-Assembly/Neuheits-Klassifikation
# ist NICHT umgesetzt.


rule rarefaction_reference_panel:
    input:
        manifest="results/panel/panel_contig_manifest.tsv",
    output:
        table="results/rarefaction/rarefaction_reference_panel.tsv",
        saturation="results/rarefaction/saturation_summary.tsv",
    params:
        genome_ids=genome_ids,
    log:
        "logs/rarefaction/reference_panel.log",
    conda:
        "../../envs/core.yaml"
    shell:
        "mkdir -p results/rarefaction logs/rarefaction && "
        "python3 workflow/scripts/rarefaction_reference_panel.py "
        "--manifest {input.manifest} "
        "--genome-ids {params.genome_ids} "
        "--out {output.table} "
        "--out-saturation {output.saturation} "
        "> {log} 2>&1"


rule extract_unmapped_reads:
    # Abschnitt 14.3, Schritt 1 (Teilumsetzung): Reads, die dem Panel
    # ueberhaupt nicht zugeordnet werden konnten.
    threads: config["threads_default"]
    input:
        bam="results/mapping/{sample_id}.panel.bam",
        bai="results/mapping/{sample_id}.panel.bam.bai",
    output:
        fastq="results/rarefaction/{sample_id}_unmapped.fastq.gz",
    conda:
        "../../envs/core.yaml"
    shell:
        "mkdir -p results/rarefaction && "
        "samtools view -@ {threads} -b -f 4 {input.bam} "
        "| samtools fastq -@ {threads} - 2>/dev/null | gzip > {output.fastq}"


rule unmapped_read_stats:
    input:
        fastq="results/rarefaction/{sample_id}_unmapped.fastq.gz",
    output:
        stats="results/rarefaction/{sample_id}_unmapped_stats.tsv",
    conda:
        "../../envs/core.yaml"
    shell:
        "seqkit stats -a {input.fastq} > {output.stats}"


rule rarefaction_all:
    input:
        "results/rarefaction/rarefaction_reference_panel.tsv",
        "results/rarefaction/saturation_summary.tsv",
        expand("results/rarefaction/{sample_id}_unmapped_stats.tsv", sample_id=mappable_sample_ids),
