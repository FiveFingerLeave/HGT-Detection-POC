# Phase VI (Section 11): window-based PAV analysis. Two mosdepth
# evaluation levels per isolate (Section 10.2):
# - "unique": MAPQ-filtered, primary alignments only (quantitative PAV)
# - "all": unfiltered, incl. secondary/supplementary (homology control,
#   reveals multi-mapping artifacts when a window appears "present"
#   only at this level)
#
# Present/Absent/Uncertain rules per region type from thresholds.yaml
# (pav.present_breadth_core/accessory, pav.absent_breadth).
#
# Phase VII (SV calling, Sniffles2) is not yet implemented.


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
    # MAPQ-filtered; mosdepth already excludes secondary/supplementary
    # by default (default --flag 1796).
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
    # --flag 1540 = exclude only unmapped(4)+qcfail(512)+dup(1024),
    # KEEP secondary(256)/supplementary - for the homology comparison.
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


# Phase VII (Section 12): breakpoint-based SV calling with Sniffles2,
# as an additional, independent line of evidence alongside the
# coverage-based PAV analysis above (Section 12: "PAV evidence =
# coverage breadth + mapping uniqueness + SV breakpoints + spanning reads").


rule sniffles_call:
    threads: config["threads_default"]
    input:
        bam="results/mapping/{sample_id}.panel.bam",
        bai="results/mapping/{sample_id}.panel.bam.bai",
        panel="results/panel/Mo_multiref_panel_v1.fa",
    output:
        vcf="results/pav/{sample_id}.sv.vcf.gz",
        snf="results/pav/{sample_id}.snf",
    log:
        "logs/pav/{sample_id}_sniffles.log",
    conda:
        "../../envs/longreads.yaml"
    shell:
        "mkdir -p results/pav logs/pav && "
        "sniffles --input {input.bam} --vcf {output.vcf} --snf {output.snf} "
        "--reference {input.panel} --threads {threads} > {log} 2>&1"


rule sniffles_cohort:
    threads: config["threads_default"]
    input:
        snf=expand("results/pav/{sample_id}.snf", sample_id=mappable_sample_ids),
    output:
        "results/pav/pilot_cohort.sv.vcf.gz",
    log:
        "logs/pav/sniffles_cohort.log",
    conda:
        "../../envs/longreads.yaml"
    shell:
        "sniffles --input {input.snf} --vcf {output} --threads {threads} > {log} 2>&1"
