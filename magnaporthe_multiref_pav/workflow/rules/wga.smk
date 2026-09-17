# Phase II (Section 7.2): Whole-genome Alignments (nucmer/MUMmer4) +
# SyRI-Klassifikation (syntenisch/umgeordnet/lokale Varianten) zwischen
# mehreren Ankerreferenzen (mind. je eine aus Oryza-, Triticum-,
# Eleusine-, Wildgras-assoziierten Linien - siehe references.tsv
# host_group-Spalte).
#
# Alle 10 ungerichteten Paare des 5-Genom-Panels (genome_pairs, siehe
# Snakefile) statt eines einzelnen Ankers - vermeidet genau die vom
# Dokument gewarnte Verzerrung durch einen einzelnen Oryza-Anker.


rule nucmer_align:
    threads: config["threads_default"]
    input:
        ref="data/references/{ref}.fa",
        query="data/references/{query}.fa",
    output:
        "results/wga/{ref}_vs_{query}.delta",
    params:
        prefix="results/wga/{ref}_vs_{query}",
    log:
        "logs/wga/{ref}_vs_{query}_nucmer.log",
    conda:
        "../../envs/wga.yaml"
    shell:
        "mkdir -p results/wga && "
        "nucmer --maxmatch -l 100 -c 500 -t {threads} "
        "-p {params.prefix} {input.ref} {input.query} > {log} 2>&1"


rule delta_filter:
    input:
        "results/wga/{ref}_vs_{query}.delta",
    output:
        "results/wga/{ref}_vs_{query}.filtered.delta",
    log:
        "logs/wga/{ref}_vs_{query}_deltafilter.log",
    conda:
        "../../envs/wga.yaml"
    shell:
        "delta-filter -m -i 85 -l 500 {input} > {output} 2> {log}"


rule show_coords:
    input:
        "results/wga/{ref}_vs_{query}.filtered.delta",
    output:
        "results/wga/{ref}_vs_{query}.coords.tsv",
    conda:
        "../../envs/wga.yaml"
    shell:
        "show-coords -THrd {input} > {output}"


rule run_syri:
    # The document's example omits -d (.delta); the installed SyRI version
    # needs it explicitly for SNP/indel identification, since the table
    # coords (-F T, default) contain no CIGAR ("CIGAR string or .delta file
    # is required"). Also, --prefix expects only the filename suffix, not a
    # path - the working directory is set separately via --dir (otherwise
    # a crash warning/error occurs).
    input:
        coords="results/wga/{ref}_vs_{query}.coords.tsv",
        delta="results/wga/{ref}_vs_{query}.filtered.delta",
        ref_fasta="data/references/{ref}.fa",
        query_fasta="data/references/{query}.fa",
    output:
        "results/syri/{ref}_vs_{query}_syri.out",
    params:
        dir="results/syri",
        prefix="{ref}_vs_{query}_",
    log:
        "logs/wga/{ref}_vs_{query}_syri.log",
    conda:
        "../../envs/wga.yaml"
    shell:
        "mkdir -p {params.dir} && "
        "syri -c {input.coords} -d {input.delta} "
        "-r {input.ref_fasta} -q {input.query_fasta} "
        "--dir {params.dir} --prefix {params.prefix} > {log} 2>&1"


rule wga_all:
    input:
        # nucmer/coords fuer alle 10 Paare, SyRI-Klassifikation nur fuer
        # die Paare mit gleicher Chromosomenzahl (syri_pairs, Snakefile).
        [f"results/wga/{ref}_vs_{query}.coords.tsv" for ref, query in genome_pairs],
        [f"results/syri/{ref}_vs_{query}_syri.out" for ref, query in syri_pairs],
