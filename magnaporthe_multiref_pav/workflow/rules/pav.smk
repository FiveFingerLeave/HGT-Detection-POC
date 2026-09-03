# Phase VI (Section 11): fensterbasierte PAV-Analyse. Zwei mosdepth-
# Auswertungsebenen pro Isolat (Abschnitt 10.2):
# - "unique": MAPQ-gefiltert, nur primaere Alignments (quantitative PAV)
# - "all": ungefiltert, inkl. secondary/supplementary (Homologie-Kontrolle,
#   deckt Multi-Mapping-Artefakte auf, wenn ein Fenster nur in dieser
#   Ebene "praesent" erscheint)
#
# Present/Absent/Uncertain-Regeln je Regionstyp aus thresholds.yaml
# (pav.present_breadth_core/accessory, pav.absent_breadth).
#
# Phase VII (SV-Calling, Sniffles2) ist noch NICHT implementiert.


rule panel_windows:
    input:
        fai="results/panel/Mo_multiref_panel_v1.fa.fai",
    output:
        genome="results/panel/panel.genome",
        windows="results/panel/panel_10kb_windows.bed",
    params:
        window_size=thresholds["pav"]["window_size_bp"],
    conda:
        "../../envs/core.yaml"
    shell:
        "cut -f1,2 {input.fai} > {output.genome} && "
        "bedtools makewindows -g {output.genome} -w {params.window_size} "
        "| awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,$3,\"win_\"NR}}' > {output.windows}"


rule mosdepth_unique:
    # MAPQ-gefiltert, mosdepth schliesst secondary/supplementary standardmaessig
    # bereits aus (Default --flag 1796).
    threads: config["threads_default"]
    input:
        bam="results/mapping/{sample_id}.panel.bam",
        bai="results/mapping/{sample_id}.panel.bam.bai",
        windows="results/panel/panel_10kb_windows.bed",
    output:
        regions="results/mapping/{sample_id}.unique.regions.bed.gz",
        thresholds_out="results/mapping/{sample_id}.unique.thresholds.bed.gz",
    params:
        prefix="results/mapping/{sample_id}.unique",
        mapq=thresholds["mapping"]["min_mapq_unique"],
        depth=thresholds["pav"]["min_depth_present"],
    log:
        "logs/pav/{sample_id}_mosdepth_unique.log",
    conda:
        "../../envs/core.yaml"
    shell:
        "mkdir -p logs/pav && "
        "mosdepth --by {input.windows} -Q {params.mapq} "
        "--thresholds {params.depth} -t {threads} "
        "{params.prefix} {input.bam} > {log} 2>&1"


rule mosdepth_all:
    # --flag 1540 = nur unmapped(4)+qcfail(512)+dup(1024) ausschliessen,
    # secondary(256)/supplementary BEHALTEN - fuer den Homologie-Vergleich.
    threads: config["threads_default"]
    input:
        bam="results/mapping/{sample_id}.panel.bam",
        bai="results/mapping/{sample_id}.panel.bam.bai",
        windows="results/panel/panel_10kb_windows.bed",
    output:
        regions="results/mapping/{sample_id}.all.regions.bed.gz",
        thresholds_out="results/mapping/{sample_id}.all.thresholds.bed.gz",
    params:
        prefix="results/mapping/{sample_id}.all",
        depth=thresholds["pav"]["min_depth_present"],
    log:
        "logs/pav/{sample_id}_mosdepth_all.log",
    conda:
        "../../envs/core.yaml"
    shell:
        "mkdir -p logs/pav && "
        "mosdepth --by {input.windows} --flag 1540 "
        "--thresholds {params.depth} -t {threads} "
        "{params.prefix} {input.bam} > {log} 2>&1"


rule pav_call:
    input:
        unique_thresholds="results/mapping/{sample_id}.unique.thresholds.bed.gz",
        unique_regions="results/mapping/{sample_id}.unique.regions.bed.gz",
        all_thresholds="results/mapping/{sample_id}.all.thresholds.bed.gz",
    output:
        windows="results/pav/{sample_id}_window_calls.tsv",
        regions="results/pav/{sample_id}_region_calls.tsv",
    params:
        present_core=thresholds["pav"]["present_breadth_core"],
        present_accessory=thresholds["pav"]["present_breadth_accessory"],
        absent=thresholds["pav"]["absent_breadth"],
    conda:
        "../../envs/core.yaml"
    shell:
        "mkdir -p results/pav && "
        "python3 workflow/scripts/pav_call.py "
        "--sample-id {wildcards.sample_id} "
        "--unique-thresholds {input.unique_thresholds} "
        "--unique-regions {input.unique_regions} "
        "--all-thresholds {input.all_thresholds} "
        "--present-breadth-core {params.present_core} "
        "--present-breadth-accessory {params.present_accessory} "
        "--absent-breadth {params.absent} "
        "--out-windows {output.windows} "
        "--out-regions {output.regions}"


rule pav_matrix:
    input:
        expand("results/pav/{sample_id}_region_calls.tsv", sample_id=mappable_sample_ids),
    output:
        "results/pav/pav_matrix.tsv",
    conda:
        "../../envs/core.yaml"
    shell:
        "python3 workflow/scripts/build_pav_matrix.py --region-calls {input} --out {output}"
