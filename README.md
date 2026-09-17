# POC: Reliable Detection of HGT Candidates and Starships in *Magnaporthe oryzae*

Proof-of-concept repository for the dissertation exposé *"Population Genomics
of Horizontal Gene Transfer in Magnaporthe oryzae"*. Authoritative workflow
specification (8 phases, assembly-/pangenome-/phylogeny-based):
[Documentation/POC_HGT_Starships_Workflow.md](Documentation/POC_HGT_Starships_Workflow.md).

This supersedes the earlier short-read-coverage-only POC approach
(see `docs/decisions.md` for the pivot rationale and what was carried over).

**Quickstart:** the actively developed, runnable part of this POC lives in
[`magnaporthe_multiref_pav/`](magnaporthe_multiref_pav/README.md), which has
its own Quickstart (environment setup, input-data population, how to run
Snakemake) and a detailed status report. License: [MIT](LICENSE).

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

Two pivots since the initial structural reset (2026-09-01). The original
8-phase `poc_hgt_starships/` skeleton below (per-phase READMEs, no
implementation) still frames the overall POC, but the actual implementation
work is happening in a dedicated sub-project,
[`magnaporthe_multiref_pav/`](magnaporthe_multiref_pav/README.md), which
elaborates the reference-panel/PAV portion of phases 1, 3, and 4 in much
more concrete detail (its own 19-section workflow document and Snakemake
pipeline).

As of 2026-09-03, all nine core sections of that sub-workflow (QC,
repeat-masking, Liftoff gene annotation, OrthoFinder orthogroups,
whole-genome alignment, Starship/Captain-gene calling, region-type
classification, panel construction, long-read mapping, PAV analysis, SV
calling, rarefaction, candidate-region manifest) have been run on a
5-genome host-representative panel against 5 real long-read test isolates
— see [`magnaporthe_multiref_pav/README.md`](magnaporthe_multiref_pav/README.md)
for the full status and headline findings (two `starship_like` regions
confirmed present across host lineages; the panel is not yet saturated at
5 genomes). BRAKER3 de-novo annotation remains blocked by a GeneMark
license, and final report synthesis (`report.smk`) is not yet implemented.

Phases 5–8 of the original 8-phase plan (validation, phylogenetic
incongruence, GLMM sample-size simulation, synthesis) are not yet started.

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
  `envs/python.yaml` remain useful tool building blocks across phases.
  Phase-specific environments (BUSCO, RepeatModeler/Masker, OrthoFinder,
  Liftoff, etc.) already exist under
  [`magnaporthe_multiref_pav/envs/`](magnaporthe_multiref_pav/envs/), with
  versions pinned to what was actually validated there — see that
  sub-project's [Quickstart](magnaporthe_multiref_pav/README.md#quickstart)
  before adding equivalents under `poc_hgt_starships/`.

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
├── magnaporthe_multiref_pav/  # actively developed panel/PAV sub-project
│                               (own Snakemake workflow, config, envs — see
│                               its README for Quickstart and status)
├── data/                   # NCBI catalogs, isolates_poc/, references/
├── assemblies/             # downloaded reference genomes (local cache)
├── envs/                   # conda environment.yml per tool group
├── docs/                   # decisions.md (chronological project log),
│                             ai_usage.md (AI usage log)
├── Documentation/          # POC_HGT_Starships_Workflow.md + planning notes
└── LICENSE                 # MIT
```

## Next step

Expand the `magnaporthe_multiref_pav/` reference panel beyond the current
5 host representatives — rarefaction analysis shows it is not yet
saturated (~19–20% candidate gain from the 5th genome) — then proceed to
Phases 5–8 of the original plan (validation, phylogenetic incongruence,
GLMM simulation, synthesis).
