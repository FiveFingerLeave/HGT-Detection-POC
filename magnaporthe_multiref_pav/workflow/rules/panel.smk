# Phase III (Section 8): Aufbau des analytischen Multi-Referenzpanels
# (Deduplikation strict-core-Bloecke, Panel-FASTA-IDs, panel_contig_manifest.tsv).
#
# Voraussetzung: panel_regions.bed aus starships.smk (Section 7.3)
# liegt vor. Der volle 14-Genom-Katalog (Panel-Ebene A, Abschnitt 8.1)
# ist bereits als config/references_full_catalog_14genomes.tsv
# archiviert; hier geht es um Panel-Ebene B (das deduplizierte
# analytische Panel aus den 5 Host-Repraesentanten).


rule build_panel:
    input:
        regions_bed="results/panel/panel_regions.bed",
        fastas=expand("data/references/{genome_id}.fa", genome_id=genome_ids),
        references_tsv=config["references_tsv"],
    output:
        manifest="results/panel/panel_contig_manifest.tsv",
        fasta="results/panel/Mo_multiref_panel_v1.fa",
    params:
        genome_ids=genome_ids,
        fasta_template=lambda wc: "data/references/{genome_id}.fa",
        min_seq_id=thresholds["panel"]["dedup_identity"],
        min_coverage=thresholds["panel"]["dedup_coverage"],
        workdir="results/panel/build_tmp",
    log:
        "logs/panel/build_panel.log",
    conda:
        "../../envs/panel.yaml"
    shell:
        "python3 workflow/scripts/build_panel.py "
        "--regions-bed {input.regions_bed} "
        "--genome-ids {params.genome_ids} "
        "--fasta-template '{params.fasta_template}' "
        "--references-tsv {input.references_tsv} "
        "--min-seq-id {params.min_seq_id} "
        "--min-coverage {params.min_coverage} "
        "--workdir {params.workdir} "
        "--out-manifest {output.manifest} "
        "--out-fasta {output.fasta} "
        "> {log} 2>&1"


rule index_panel_fasta:
    input:
        "results/panel/Mo_multiref_panel_v1.fa",
    output:
        "results/panel/Mo_multiref_panel_v1.fa.fai",
    conda:
        "../../envs/panel.yaml"
    shell:
        "samtools faidx {input}"
