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


# Phase VIII (Section 13): Kandidatenregionen mit stabilen IDs versehen
# (13.1) und pro Testisolat die kombinierte PAV+SV-Evidenz zuordnen
# (13.2). Abschnitt 13.3 (bedtools-Intersect mit Gen-GFF3 fuer
# Funktionsannotation je PAV-Block) ist NICHT umgesetzt - dafuer
# fehlt eine einheitliche Panel-weite Gen-GFF3 (unsere Genannotation
# liegt nur pro Ausgangsgenom vor, nicht auf Panel-Koordinaten
# projiziert).


rule build_candidate_regions:
    input:
        manifest="results/panel/panel_contig_manifest.tsv",
        yr_gff="results/starships/panel_YR.filt.gff",
        starship_candidates="results/starships/starship_like_candidates.tsv",
    output:
        "results/panel/panel_candidate_regions.tsv",
    log:
        "logs/panel/build_candidate_regions.log",
    conda:
        "../../envs/core.yaml"
    shell:
        "python3 workflow/scripts/build_candidate_regions.py "
        "--manifest {input.manifest} "
        "--yr-gff {input.yr_gff} "
        "--starship-candidates {input.starship_candidates} "
        "--out {output} > {log} 2>&1"


rule build_candidate_region_calls:
    input:
        candidates="results/panel/panel_candidate_regions.tsv",
        window_calls=expand("results/pav/{sample_id}_window_calls.tsv", sample_id=mappable_sample_ids),
        region_calls=expand("results/pav/{sample_id}_region_calls.tsv", sample_id=mappable_sample_ids),
        sv_vcf="results/pav/pilot_cohort.sv.vcf.gz",
    output:
        "results/pav/candidate_region_calls.tsv",
    log:
        "logs/panel/build_candidate_region_calls.log",
    conda:
        "../../envs/core.yaml"
    shell:
        "python3 workflow/scripts/build_candidate_region_calls.py "
        "--candidates {input.candidates} "
        "--window-calls {input.window_calls} "
        "--region-calls {input.region_calls} "
        "--sv-vcf {input.sv_vcf} "
        "--out {output} > {log} 2>&1"
