# Section 3 of the POC guideline: core vs. candidate PAV matrices, PCA +
# k-means clustering, Adjusted Rand Index against host-lineage labels.

rule combine_host_labels:
    # Short-read (config/samples.tsv) and long-read (config/samples_longread.tsv)
    # isolates carry host_lineage in separate tables; clustering needs one
    # unified sample_id -> host_lineage lookup covering every mapped sample.
    input:
        short="config/samples.tsv",
        long="config/samples_longread.tsv",
    output:
        "results/clustering/host_labels.tsv",
    run:
        import pandas as pd

        short_df = pd.read_csv(input.short, sep="\t", dtype=str)[["sample_id", "host_lineage"]]
        long_df = pd.read_csv(input.long, sep="\t", dtype=str)[["sample_id", "host_lineage"]]
        pd.concat([short_df, long_df], ignore_index=True).to_csv(output[0], sep="\t", index=False)


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
        host_labels="results/clustering/host_labels.tsv",
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
