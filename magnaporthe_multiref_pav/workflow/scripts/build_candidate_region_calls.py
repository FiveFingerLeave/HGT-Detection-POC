#!/usr/bin/env python3
"""Section 13.2: per-isolate PAV status for each candidate region,
combining the window-level PAV calls (Section 11) with SV support
(Section 12, Sniffles2 cohort VCF)."""
import argparse
import csv
import gzip
from collections import defaultdict


def load_candidates(path):
    with open(path) as f:
        return {r["candidate_id"]: r for r in csv.DictReader(f, delimiter="\t")}


def load_window_calls(paths):
    """sample_id -> panel_id -> list of (breadth_unique, mean_depth, is_ambiguous)"""
    data = defaultdict(lambda: defaultdict(list))
    for path in paths:
        with open(path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                data[r["sample_id"]][r["panel_id"]].append((
                    float(r["breadth_unique"]), float(r["mean_depth_unique"]),
                    r["call"] == "ambiguous_multimapping",
                ))
    return data


def load_region_calls(paths):
    """sample_id -> panel_id -> call"""
    data = defaultdict(dict)
    for path in paths:
        with open(path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                data[r["sample_id"]][r["panel_id"]] = r["call"]
    return data


def load_sv_support(vcf_path, candidate_chroms):
    """chrom -> sample_id -> count of non-ref genotype calls."""
    support = defaultdict(lambda: defaultdict(int))
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        samples = []
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            cols = line.rstrip("\n").split("\t")
            chrom = cols[0]
            if chrom not in candidate_chroms:
                continue
            for sample, gt_field in zip(samples, cols[9:]):
                gt = gt_field.split(":")[0]
                if gt not in ("./.", "0/0", "."):
                    support[chrom][sample] += 1
    return support


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--window-calls", nargs="+", required=True)
    ap.add_argument("--region-calls", nargs="+", required=True)
    ap.add_argument("--sv-vcf", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    candidates = load_candidates(args.candidates)
    window_data = load_window_calls(args.window_calls)
    region_calls = load_region_calls(args.region_calls)
    candidate_chroms = {c["panel_id"] for c in candidates.values()}

    # Window/region calls only store panel_id, but the SV VCF's CHROM is
    # the full panel FASTA header - recover the mapping from the VCF body.
    chrom_lookup = {}
    opener = gzip.open if args.sv_vcf.endswith(".gz") else open
    with opener(args.sv_vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0]
            panel_id = chrom.split("|", 1)[0]
            if panel_id in candidate_chroms and panel_id not in chrom_lookup:
                chrom_lookup[panel_id] = chrom
            if len(chrom_lookup) == len(candidate_chroms):
                break

    full_chroms = set(chrom_lookup.values())
    sv_support = load_sv_support(args.sv_vcf, full_chroms)

    samples = sorted(region_calls.keys())
    rows = []
    for candidate_id, cand in candidates.items():
        panel_id = cand["panel_id"]
        full_chrom = chrom_lookup.get(panel_id)
        for sample in samples:
            status = region_calls[sample].get(panel_id, "NA")
            windows = window_data[sample].get(panel_id, [])
            n_windows = len(windows)
            breadth_5x = sum(w[0] for w in windows) / n_windows if n_windows else 0.0
            mean_depth = sum(w[1] for w in windows) / n_windows if n_windows else 0.0
            n_ambiguous = sum(1 for w in windows if w[2])
            unique_window_fraction = 1 - (n_ambiguous / n_windows) if n_windows else 0.0
            sv_n = sv_support.get(full_chrom, {}).get(sample, 0) if full_chrom else 0

            if status == "present" and unique_window_fraction >= 0.8:
                confidence = "high" if sv_n > 0 else "medium"
            elif status == "present":
                confidence = "low"
            elif status == "absent":
                confidence = "high"
            else:
                confidence = "low"

            notes = []
            if status == "ambiguous_multimapping":
                notes.append("repeat-associated multi-mapping")
            if cand["region_type"] == "starship_like" and cand["captain_gene_id"] != "NA":
                notes.append(f"captain={cand['captain_gene_id']}")

            rows.append({
                "test_isolate": sample,
                "candidate_id": candidate_id,
                "region_type": cand["region_type"],
                "status": status,
                "confidence": confidence,
                "breadth_5x": f"{breadth_5x:.4f}",
                "mean_depth": f"{mean_depth:.2f}",
                "unique_window_fraction": f"{unique_window_fraction:.4f}",
                "sv_support": sv_n,
                "shared_reference_isolates": cand["reference_isolates"],
                "closest_panel_reference": cand["representative_reference"],
                "assignment_confidence": confidence,
                "notes": ";".join(notes) if notes else "NA",
            })

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} candidate-region x isolate calls to {args.out}")


if __name__ == "__main__":
    main()
