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
        ari_summary="results/clustering/ari_summary.tsv",
        # Per-region ARI (candidate matrix only): a more sensitive,
        # targeted complement to the whole-matrix ARI above, since a
        # single discordant region can be diluted out of an aggregate
        # PCA/k-means score by several unrelated candidate regions - see
        # docs/decisions.md.
        per_region_ari="results/clustering/candidate_per_region_ari.tsv",
    params:
        k=config["clustering"]["k"],
        n_components=config["clustering"]["n_components"],
        z_threshold=config["clustering"]["discordance_z_threshold"],
    log:
        "logs/clustering/cluster_and_score.log",
    conda:
        "../../envs/python.yaml"
    shell:
        "python workflow/scripts/cluster_and_score.py "
        "--core {input.core} --candidate {input.candidate} "
        "--host-labels {input.host_labels} "
        "--k {params.k} --n-components {params.n_components} "
        "--discordance-z-threshold {params.z_threshold} "
        "--output {output.ari_summary} "
        "--per-region-output {output.per_region_ari} "
        "> {log} 2>&1"
