# POC: Horizontal Gene Transfer in *Magnaporthe oryzae*

Proof-of-concept repository for the dissertation exposé *"Population Genomics
of Horizontal Gene Transfer in Magnaporthe oryzae"* (WP1.2–WP1.4). Full
methodological guideline:
[Dokumentation/POC Workflow – HGT Machbarkeitsstudie M. oryzae.md](Dokumentation/POC%20Workflow%20%E2%80%93%20HGT%20Machbarkeitsstudie%20M.%20oryzae.md).

## POC scope (from the guideline)

1. **Multireference panel functionality** — short reads mapped against a
   4–6 genome panel (host lineages + *M. grisea* outgroup), coverage/PAV
   candidate calling.
2. **Starship/HGT candidate counting** — coverage-based (Step 1) +
   structure-based (Starfish/Stargraph), cross-validated.
3. **Core vs. candidate clustering by host type** — PCA/ARI discordance
   test; the diagnostic signature of HGT.
4. **GLMM power analysis (Q3)** — simulation-based, data-independent,
   runnable now.

## Current status

No real isolate short-read data has been pulled into this repo yet.
Steps 1–3 are executable protocols, ready to run once reference genomes
and an isolate subset are staged under `data/`. Step 4 can be run now
(`power_analysis/glmm_power_sim.R`) with the placeholder effect sizes from
the guideline.

## Data

- `data/ncbi_m_oryzae_sra_runs.tsv` — full NCBI SRA catalog for
  *Pyricularia oryzae* (3,902 runs; all library strategies/platforms).
- `data/ncbi_m_oryzae_sra_wgs_illumina.tsv` — WGS short-read subset
  (1,754 runs) — candidate pool for the POC isolate subset (Step 1).
- `data/ncbi_m_oryzae_sra_wgs_longread.tsv` — WGS long-read subset
  (193 runs, Nanopore/PacBio) — candidate pool for reference-panel or
  future high-quality assemblies.
- `data/ncbi_m_oryzae_assemblies.tsv`/`.jsonl` — full NCBI assembly-level
  catalog (finished/submitted genomes, not raw reads).
- `data/ncbi_m_oryzae_longread_candidates*.tsv`,
  `ncbi_m_oryzae_longread_highquality.tsv` — pre-filtered high-quality
  assembly candidates.
- `data/assembly_accessions_selected.txt`, `starfish_panel_candidates.tsv`
  — earlier accession shortlists (predate this reset; re-validate before
  reuse).
- `input/metadata/ncbi_isolates_host.tsv` — curated host species per
  assembly accession.
- `assemblies/` — previously downloaded NCBI assembly FASTA/GFF (local
  cache, not tracked in git; see `.gitignore`).
- `data/references/`, `data/isolates_poc/` — empty, to be populated with
  the actual multireference panel and POC isolate short reads.

## Repo structure

```text
barragan/
├── data/
│   ├── references/        # multireference panel (FASTA + GFF)
│   ├── isolates_poc/      # Illumina/Nanopore subset for the POC
│   └── *.tsv/.jsonl        # NCBI isolate/sequence catalogs (see above)
├── envs/                  # one conda environment.yml per tool group
├── workflow/
│   ├── Snakefile
│   ├── rules/             # mapping.smk, pav.smk, starfish.smk, cluster.smk
│   └── scripts/
├── results/
│   ├── mapping/
│   ├── pav_calls/
│   ├── starfish/
│   └── clustering/
├── power_analysis/
│   └── glmm_power_sim.R
├── docs/                  # decisions.md, ai_usage.md
└── Dokumentation/         # full guideline + prior planning notes
```

## Next steps

See "Zusammenfassung: Offene Punkte für die Umsetzung in VS Code" in the
guideline document.
