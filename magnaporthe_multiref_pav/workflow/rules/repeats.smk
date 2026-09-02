# Phase I (Section 6.3): RepeatModeler2 + RepeatMasker auf allen 14
# Referenzen; Repeat-Anteil pro Contig/Analysefenster.
#
# Repeatreiche Bereiche werden NICHT entfernt, sondern nur annotiert - die
# spaetere PAV-Klassifikation (pav.smk) markiert Fenster mit hohem
# Repeat-Anteil als moeglicherweise mehrdeutig, statt sie allein wegen
# geringer eindeutiger Coverage als "absent" zu werten.
#
# Eigene Environment (envs/repeats.yaml) statt der schweren
# envs/annotation.yaml (BRAKER3 etc.) - RepeatModeler/RepeatMasker
# brauchen keine der dortigen Abhaengigkeiten.
#
# Laufzeit-Hinweis: RepeatModeler mit -LTRStruct ist der teuerste Schritt
# in Phase I (Stunden pro Genom, siehe docs/decisions.md fuer den
# Pilotlauf-Zeitwert). Wie bei busco_reference (qc.smk) reserviert
# threads: das volle Kernbudget pro Job, damit Snakemake bei --cores 1
# echte Serialisierung erzwingt statt mehrere RepeatModeler-Laeufe
# gleichzeitig zuzulassen (gleiches OOM-Risiko wie bei BUSCO).


rule index_reference_fai:
    input:
        "data/references/{genome_id}.fa",
    output:
        "data/references/{genome_id}.fa.fai",
    conda:
        "../../envs/core.yaml"
    shell:
        "samtools faidx {input}"


rule build_repeat_database:
    input:
        fasta="data/references/{genome_id}.fa",
    output:
        "results/repeats/{genome_id}/{genome_id}.nin",
    params:
        out_dir="results/repeats/{genome_id}",
    log:
        "logs/repeats/{genome_id}_builddatabase.log",
    conda:
        "../../envs/repeats.yaml"
    shell:
        "mkdir -p {params.out_dir} && "
        "cd {params.out_dir} && "
        "BuildDatabase -name {wildcards.genome_id} -engine ncbi ../../../{input.fasta} "
        "> ../../../{log} 2>&1"


rule run_repeatmodeler:
    threads: config["threads_default"]
    input:
        "results/repeats/{genome_id}/{genome_id}.nin",
    output:
        "results/repeats/{genome_id}/{genome_id}-families.fa",
    params:
        out_dir="results/repeats/{genome_id}",
    log:
        "logs/repeats/{genome_id}_repeatmodeler.log",
    conda:
        "../../envs/repeats.yaml"
    shell:
        "cd {params.out_dir} && "
        "RepeatModeler -database {wildcards.genome_id} -threads {threads} -LTRStruct "
        "> ../../../{log} 2>&1"


rule run_repeatmasker:
    threads: config["threads_default"]
    input:
        fasta="data/references/{genome_id}.fa",
        lib="results/repeats/{genome_id}/{genome_id}-families.fa",
    output:
        out="results/repeats/{genome_id}/masked/{genome_id}.fa.out",
        masked="results/repeats/{genome_id}/masked/{genome_id}.fa.masked",
        gff="results/repeats/{genome_id}/masked/{genome_id}.fa.out.gff",
    params:
        out_dir="results/repeats/{genome_id}/masked",
    log:
        "logs/repeats/{genome_id}_repeatmasker.log",
    conda:
        "../../envs/repeats.yaml"
    shell:
        "mkdir -p {params.out_dir} && "
        "RepeatMasker -pa {threads} -lib {input.lib} -gff -dir {params.out_dir} "
        "{input.fasta} > {log} 2>&1"


rule repeat_windows:
    input:
        out="results/repeats/{genome_id}/masked/{genome_id}.fa.out",
        fai="data/references/{genome_id}.fa.fai",
    output:
        windows="results/repeats/{genome_id}_repeat_windows.bed",
        per_contig="results/repeats/{genome_id}_repeat_per_contig.tsv",
    params:
        window_size=thresholds["pav"]["window_size_bp"],
    conda:
        "../../envs/repeats.yaml"
    shell:
        "bash workflow/scripts/repeat_windows.sh "
        "{input.out} {input.fai} {params.window_size} "
        "{output.windows} {output.per_contig}"


rule repeats_all:
    input:
        expand("results/repeats/{genome_id}_repeat_windows.bed", genome_id=genome_ids),
        expand("results/repeats/{genome_id}_repeat_per_contig.tsv", genome_id=genome_ids),
