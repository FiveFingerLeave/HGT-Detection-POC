# Section 3 of the POC guideline: core vs. candidate PAV matrices, PCA +
# k-means clustering, Adjusted Rand Index against host-lineage labels.

rule build_pav_matrices:
    input:
        candidate_table="results/pav_calls/candidate_table.tsv",
    output:
        core="results/clustering/core_pav_matrix.tsv",
        candidate="results/clustering/candidate_pav_matrix.tsv",
    log:
        "logs/clustering/build_pav_matrices.log",
    conda:
        "../../envs/python.yaml"
    shell:
        "python workflow/scripts/build_pav_matrix.py "
        "--candidate-table {input.candidate_table} "
        "--output-core {output.core} "
        "--output-candidate {output.candidate} "
        "> {log} 2>&1"


rule cluster_and_score:
    input:
        core="results/clustering/core_pav_matrix.tsv",
        candidate="results/clustering/candidate_pav_matrix.tsv",
        host_labels="config/samples.tsv",
    output:
        "results/clustering/ari_summary.tsv",
    params:
        k=config["clustering"]["k"],
        n_components=config["clustering"]["n_components"],
    log:
        "logs/clustering/cluster_and_score.log",
    conda:
        "../../envs/python.yaml"
    shell:
        "python workflow/scripts/cluster_and_score.py "
        "--core {input.core} --candidate {input.candidate} "
        "--host-labels {input.host_labels} "
        "--k {params.k} --n-components {params.n_components} "
        "--output {output} "
        "> {log} 2>&1"
