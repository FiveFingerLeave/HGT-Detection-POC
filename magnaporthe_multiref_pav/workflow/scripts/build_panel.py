#!/usr/bin/env python3
"""Section 8: build the deduplicated analytical multi-reference panel.

Pipeline:
1. Merge adjacent same-class 10kb windows (panel_regions.bed, Section
   7.3) into contiguous per-genome region blocks.
2. Keep only panel-relevant classes (8.1): strict_core, soft_core,
   shell, private_accessory, accessory_chromosome, starship_like.
   repeat_ambiguous/unclassified/subtelomeric_dynamic are excluded -
   too uncertain to serve as panel anchors.
3. Extract each block's sequence (samtools faidx, batched per genome).
4. Cluster ALL block sequences together with mmseqs2 easy-cluster at
   98% identity / 90% mutual coverage (8.2) - this single step
   implements BOTH halves of the document's dedup rule at once:
   near-identical copies of the same strict_core block across genomes
   collapse to one representative, while genuinely distinct
   soft_core/shell/accessory/starship_like variants remain separate
   clusters (no special-casing needed per region type).
5. Write the panel manifest (8.4) and the final panel FASTA with
   PANEL<n> headers (8.3), using mmseqs2's own representative choice.
"""
import argparse
import csv
import subprocess
from collections import defaultdict, Counter


PANEL_CLASSES = {
    "strict_core", "soft_core", "shell", "private_accessory",
    "accessory_chromosome", "starship_like",
}

CLUSTER_PREFIX = {
    "strict_core": "CORE", "soft_core": "CORE", "shell": "CORE",
    "private_accessory": "ACC", "accessory_chromosome": "ACC",
    "starship_like": "STAR",
}


def merge_blocks(regions_bed):
    rows = []
    with open(regions_bed) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rows.append(r)

    blocks = []
    current = None
    for r in rows:
        g, c, s, e, cls = r["genome_id"], r["contig"], int(r["start"]), int(r["end"]), r["region_class"]
        if (current and current["genome_id"] == g and current["contig"] == c
                and current["region_class"] == cls and current["end"] == s):
            current["end"] = e
            current["repeat_fracs"].append(float(r["repeat_fraction"]))
        else:
            if current:
                blocks.append(current)
            current = {"genome_id": g, "contig": c, "start": s, "end": e,
                       "region_class": cls, "repeat_fracs": [float(r["repeat_fraction"])]}
    if current:
        blocks.append(current)
    return [b for b in blocks if b["region_class"] in PANEL_CLASSES]


def block_id(b):
    return f"{b['contig']}:{b['start'] + 1}-{b['end']}"


def extract_sequences(blocks, genome_ids, fasta_template, workdir):
    combined_fasta = f"{workdir}/all_blocks.fa"
    with open(combined_fasta, "w") as out:
        for g in genome_ids:
            g_blocks = [b for b in blocks if b["genome_id"] == g]
            if not g_blocks:
                continue
            regions_file = f"{workdir}/{g}_regions.txt"
            with open(regions_file, "w") as rf:
                for b in g_blocks:
                    rf.write(block_id(b) + "\n")
            result = subprocess.run(
                ["samtools", "faidx", fasta_template.format(genome_id=g), "-r", regions_file],
                capture_output=True, text=True, check=True,
            )
            out.write(result.stdout)
    return combined_fasta


def run_mmseqs_cluster(combined_fasta, workdir, min_seq_id, cov):
    prefix = f"{workdir}/panel_clu"
    tmp = f"{workdir}/mmseqs_tmp"
    subprocess.run(
        ["mmseqs", "easy-cluster", combined_fasta, prefix, tmp,
         "--min-seq-id", str(min_seq_id), "-c", str(cov), "--cov-mode", "0"],
        check=True,
    )
    cluster_tsv = f"{prefix}_cluster.tsv"
    rep_fasta = f"{prefix}_rep_seq.fasta"
    clusters = defaultdict(list)
    with open(cluster_tsv) as f:
        for line in f:
            rep, member = line.rstrip("\n").split("\t")
            clusters[rep].append(member)
    rep_seqs = {}
    with open(rep_fasta) as f:
        name, seq = None, []
        for line in f:
            if line.startswith(">"):
                if name:
                    rep_seqs[name] = "".join(seq)
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name:
            rep_seqs[name] = "".join(seq)
    return clusters, rep_seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions-bed", required=True)
    ap.add_argument("--genome-ids", nargs="+", required=True)
    ap.add_argument("--fasta-template", required=True)
    ap.add_argument("--references-tsv", required=True)
    ap.add_argument("--min-seq-id", type=float, required=True)
    ap.add_argument("--min-coverage", type=float, required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out-manifest", required=True)
    ap.add_argument("--out-fasta", required=True)
    args = ap.parse_args()

    import os
    os.makedirs(args.workdir, exist_ok=True)

    host_group = {}
    with open(args.references_tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            host_group[r["genome_id"]] = r["host_group"]

    blocks = merge_blocks(args.regions_bed)
    # block_id() already returns "{genome_id}__{contig}:{start}-{end}" since
    # b['contig'] itself carries the genome prefix (our fasta header
    # convention) - this is also exactly the sequence name samtools faidx
    # and mmseqs2 report, so it's used directly as the lookup key.
    block_lookup = {block_id(b): b for b in blocks}

    combined_fasta = extract_sequences(blocks, args.genome_ids, args.fasta_template, args.workdir)
    clusters, rep_seqs = run_mmseqs_cluster(combined_fasta, args.workdir, args.min_seq_id, args.min_coverage)

    def lookup_block(member_header):
        return block_lookup[member_header]

    cluster_counters = Counter()
    manifest_rows = []
    panel_records = []

    cluster_items = []
    for rep, members in clusters.items():
        member_blocks = [lookup_block(m) for m in members]
        class_counts = Counter(b["region_class"] for b in member_blocks)
        region_type = class_counts.most_common(1)[0][0]
        rep_block = lookup_block(rep)
        cluster_items.append((region_type, rep_block["genome_id"], rep_block["start"], rep, members, member_blocks))

    cluster_items.sort(key=lambda x: (x[0] != "strict_core", x[0], x[1], x[2]))

    for idx, (region_type, _, _, rep, members, member_blocks) in enumerate(cluster_items, start=1):
        panel_id = f"PANEL{idx:06d}"
        prefix = CLUSTER_PREFIX[region_type]
        cluster_counters[prefix] += 1
        region_cluster = f"{prefix}_{cluster_counters[prefix]:03d}"

        rep_genome, rep_source = rep.split("__", 1)
        rep_block = lookup_block(rep)
        seq = rep_seqs[rep]

        genomes_in_cluster = sorted({b["genome_id"] for b in member_blocks})
        hosts_in_cluster = sorted({host_group.get(g, "unknown") for g in genomes_in_cluster})
        mean_repeat = sum(sum(b["repeat_fracs"]) / len(b["repeat_fracs"]) for b in member_blocks) / len(member_blocks)

        manifest_rows.append({
            "panel_id": panel_id,
            "region_cluster": region_cluster,
            "region_type": region_type,
            "representative_reference": rep_genome,
            "source_interval": f"{rep_genome}__{rep_source}",
            "length_bp": len(seq),
            "n_reference_genomes": len(genomes_in_cluster),
            "reference_isolates": ";".join(genomes_in_cluster),
            "host_groups": ";".join(hosts_in_cluster),
            "repeat_fraction": f"{mean_repeat:.4f}",
            "starship_evidence": "yes" if region_type == "starship_like" else "no",
        })
        panel_records.append((panel_id, region_type, region_cluster, rep_genome, rep_source, seq))

    with open(args.out_manifest, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    with open(args.out_fasta, "w") as f:
        for panel_id, region_type, region_cluster, rep_genome, rep_source, seq in panel_records:
            header = f">{panel_id}|type={region_type}|cluster={region_cluster}|rep={rep_genome}|source={rep_genome}__{rep_source}"
            f.write(header + "\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + "\n")

    print(f"Panel built: {len(panel_records)} regions from {len(blocks)} pre-cluster blocks")
    type_counts = Counter(r["region_type"] for r in manifest_rows)
    for k, v in type_counts.most_common():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
