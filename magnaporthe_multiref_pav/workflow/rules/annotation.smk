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


# Section 7.1: orthogroups. gffread extracts the protein sequences from
# the Liftoff GFF3 files (a prerequisite for orthofinder -f), then
# orthofinder itself exactly per the document's command (thread counts
# adapted to the 14 available cores instead of the document's example
# values of 32/16).


rule extract_proteome:
    # gffread -y marks stop codons/codons that can't be cleanly translated
    # (e.g. at exon boundaries with a frame remainder) with "." - diamond
    # does not accept that in its sequence alphabet (see docs/decisions.md:
    # diamond v2.2.6 even hangs completely on such a character instead of
    # raising an error). Replace "." (and, as a precaution, "-") in
    # sequence lines (not headers) with "X" (unknown amino acid).
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
        # The document specifies "-M msa -T iqtree". Tested and then
        # discarded (see docs/decisions.md): with only 5 species, the
        # orthogroup assignment itself (diamond+MCL) stays fast (minutes),
        # but the subsequent gene-tree inference runs iqtree3 WITH
        # ModelFinder PRO per orthogroup individually - with >11,000
        # orthogroups and >15 min per tree, that would be a multi-day to
        # multi-week run. For the output actually needed in Section 7.1
        # (the Orthogroups.tsv presence/absence matrix for the gene-based
        # PAV classification), that is not necessary - "-M dendroblast"
        # delivers the same orthogroup matrix (determined identically by
        # diamond+MCL) without the expensive MSA/gene-tree refinement
        # layer.
        "rm -rf {params.proteome_dir}/OrthoFinder {params.fixed_dir} && "
        "mkdir -p results/orthofinder && "
        "orthofinder -f {params.proteome_dir} -S diamond -M dendroblast "
        "-t {threads} -a {threads} > {log} 2>&1 && "
        "mv {params.proteome_dir}/OrthoFinder/Results_*/ {params.fixed_dir}"
