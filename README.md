# POC: Reliable Detection of HGT Candidates and Starships in *Magnaporthe oryzae*

Proof-of-concept repository for the dissertation exposé *"Population Genomics
of Horizontal Gene Transfer in Magnaporthe oryzae"*. Authoritative workflow
specification (8 phases, assembly-/pangenome-/phylogeny-based):
[Documentation/POC_HGT_Starships_Workflow.md](Documentation/POC_HGT_Starships_Workflow.md).

This supersedes the earlier short-read-coverage-only POC approach
(see `docs/decisions.md` for the pivot rationale and what was carried over).

## POC scope (from the guideline)

Not exact HGT frequency — proof that the workflow reliably identifies HGT/
Starship candidates, plus a data-driven recommendation for sample size and
GLMM model complexity for the dissertation's main study.

1. **Phase 1 — Reference panel**: 20–30 isolates across clone lineages/hosts,
   prioritizing existing long-read assemblies; BUSCO/N50/BlobTools QC;
   funannotate/BRAKER annotation; RepeatModeler/RepeatMasker masking.
2. **Phase 2 — Starship identification**: `starfish` per genome, manually
   curated catalog (coordinates, Captain family, cargo, insertion site).
3. **Phase 3 — Presence/absence matrix**: pangenome tool (Panaroo/PIRATE/
   Roary) + alignment-based accessory-region detection, joined with the
   Starship catalog.
4. **Phase 4 — Saturation/rarefaction analysis** (core of the POC): does
   the candidate-discovery curve flatten with panel size?
5. **Phase 5 — Validation**: leave-one-out against long-read ground truth +
   real short-read test isolates.
6. **Phase 6 — Phylogenetic incongruence**: core-genome tree vs.
   Starship-gene tree topology comparison as HGT evidence.
7. **Phase 7 — Sample-size/model-complexity simulation for the later GLMM**
   — deliberately deferred until real candidate rates exist (see
   `docs/decisions.md`).
8. **Phase 8 — Synthesis & go/no-go report.**

## Current status

Structural reset only (2026-09-01): old short-read-coverage pipeline
(`workflow/`, `config/`, old `power_analysis/`) removed; `poc_hgt_starships/`
phase skeleton created with per-phase READMEs (input/output spec, tools).
No phase is implemented yet. Phase 1 (reference-panel curation) is the
natural next step — the NCBI catalogs below already narrow the search.

## Data (kept from the previous phase, still directly useful)

- `data/ncbi_m_oryzae_sra_runs.tsv` — full NCBI SRA catalog for
  *Pyricularia oryzae* (3,902 runs; all library strategies/platforms).
- `data/ncbi_m_oryzae_sra_wgs_illumina.tsv` — WGS short-read subset
  (1,754 runs) — candidate pool for Phase 5's real short-read test isolates.
- `data/ncbi_m_oryzae_sra_wgs_longread.tsv` — WGS long-read subset
  (193 runs, Nanopore/PacBio) — candidate pool for Phase 1's reference panel.
- `data/ncbi_m_oryzae_assemblies.tsv`/`.jsonl` — full NCBI assembly-level
  catalog (finished/submitted genomes).
- `data/ncbi_m_oryzae_longread_candidates*.tsv`,
  `ncbi_m_oryzae_longread_highquality.tsv` — pre-filtered high-quality
  assembly candidates for Phase 1.
- `input/metadata/ncbi_isolates_host.tsv` — curated host species per
  assembly accession.
- `assemblies/` — previously downloaded NCBI assembly FASTA/GFF (local
  cache, not tracked in git; 6 genomes from the old panel, a starting
  subset for the new 20–30-genome Phase 1 panel).
- `data/isolates_poc/` — 5 short-read + 1 long-read isolate FASTQ already
  downloaded and QC-checked; candidates for Phase 5 test isolates (local
  cache, not tracked in git).
- `envs/starfish.yaml` — still directly reusable for Phase 2.
  `envs/mapping.yaml` (bwa-mem2/minimap2/samtools/mosdepth/fastp) and
  `envs/python.yaml` remain useful tool building blocks across phases;
  phase-specific environments (BUSCO, funannotate, RepeatMasker, Panaroo,
  IQ-TREE, etc.) still need to be added under `poc_hgt_starships/`.

## Repo structure

```text
barragan/
├── poc_hgt_starships/
│   ├── 00_data/            # panel curation, metadata
│   ├── 01_assembly_qc/     # assembly, BUSCO, BlobTools
│   ├── 02_annotation/      # funannotate/BRAKER, RepeatModeler/Masker
│   ├── 03_starship_calls/  # starfish output, curated catalog
│   ├── 04_pav_matrix/      # pangenome PAV + accessory-region detection
│   ├── 05_saturation/      # rarefaction/saturation analysis
│   ├── 06_validation/      # leave-one-out + real test isolates
│   ├── 07_phylogeny/       # core tree vs. Starship-gene tree
│   ├── 08_simulation_glmm/ # deferred - sample-size/model simulation
│   └── 09_report/          # synthesis, go/no-go
├── data/                   # NCBI catalogs, isolates_poc/, references/
├── assemblies/             # downloaded reference genomes (local cache)
├── envs/                   # conda environment.yml per tool group
├── docs/                   # decisions.md (chronological project log)
└── Documentation/          # POC_HGT_Starships_Workflow.md + planning notes
```

## Next step

Phase 1 (`poc_hgt_starships/00_data/`): curate the 20–30-isolate reference
panel from the existing NCBI long-read-assembly catalogs, maximizing clone
lineage/host diversity.
