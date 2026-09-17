# Phase III (Section 8): Building the analytical multi-reference panel
# (deduplication of strict-core blocks, panel FASTA IDs, panel_contig_manifest.tsv).
#
# Prerequisite: panel_regions.bed from starships.smk (Section 7.3)
# must be available. The full 14-genome catalog (Panel Level A, Section 8.1)
# is already archived as config/references_full_catalog_14genomes.tsv;
# this file concerns Panel Level B (the deduplicated
# analytical panel built from the 5 host representatives).


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


# Phase VIII (Section 13): Assign stable IDs to candidate regions
# (13.1) and, per test isolate, assign the combined PAV+SV evidence
# (13.2). Section 13.3 (bedtools intersect with gene GFF3 for
# functional annotation per PAV block) is NOT implemented - this is
# missing a unified panel-wide gene GFF3 (our gene annotation
# is only available per source genome, not projected onto
# panel coordinates).


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
