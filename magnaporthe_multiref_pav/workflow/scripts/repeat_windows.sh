#!/bin/bash
# Section 6.3: repeat fraction per contig and per analysis window, derived
# from a RepeatMasker .out. Uses the same window size as the
# later window-based PAV classification (config/thresholds.yaml:
# pav.window_size_bp), so the two can later be joined directly.
set -euo pipefail

RM_OUT=$1
FAI=$2
WINDOW=$3
OUT_WINDOWS=$4
OUT_CONTIG=$5

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

cut -f1,2 "$FAI" > "$TMP/genome.txt"

# RepeatMasker .out: 3 header lines, then whitespace-separated columns;
# col5=query seq, col6=begin(1-based), col7=end -> BED (0-based half-open)
awk 'NR>3 {print $5"\t"($6-1)"\t"$7}' "$RM_OUT" | sort -k1,1 -k2,2n > "$TMP/repeats_raw.bed"
bedtools merge -i "$TMP/repeats_raw.bed" > "$TMP/repeats_merged.bed"

bedtools makewindows -g "$TMP/genome.txt" -w "$WINDOW" > "$TMP/windows.bed"
bedtools coverage -a "$TMP/windows.bed" -b "$TMP/repeats_merged.bed" > "$OUT_WINDOWS"

awk 'BEGIN{FS=OFS="\t"} {print $1, ($3-$2)}' "$TMP/repeats_merged.bed" \
    | awk 'BEGIN{FS=OFS="\t"} {sum[$1]+=$2} END{for (c in sum) print c, sum[c]}' \
    > "$TMP/covered_per_contig.tsv"

awk 'BEGIN{FS=OFS="\t"}
     FNR==NR{len[$1]=$2; order[++n]=$1; next}
     {cov[$1]=$2}
     END{for (i=1; i<=n; i++) {
           c=order[i]; l=len[c]; rb=(c in cov)?cov[c]:0;
           frac=(l>0)?rb/l:0;
           print c, l, rb, frac
         }}' "$TMP/genome.txt" "$TMP/covered_per_contig.tsv" > "$TMP/body.tsv"

{ echo -e "contig\tlength_bp\trepeat_bp\trepeat_fraction"; sort -k1,1 "$TMP/body.tsv"; } > "$OUT_CONTIG"
