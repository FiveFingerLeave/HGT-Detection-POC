# Multireference Panel & Long-Read PAV Workflow

Implementation of
[Documentation/multireference_panel_pav_workflow.md](../Documentation/multireference_panel_pav_workflow.md)
(19 sections: QC, annotation, pangenome classification, panel construction,
long-read mapping, window-based PAV, SV calling, rarefaction).

## Status (2026-09-02)

**Phase I (Section 6.1, QC) — complete for the full 14-genome catalog:**
- 14 Complete-Genome references identified from the full NCBI catalog
  (`../poc_hgt_starships/00_data/ncbi_pyricularia_oryzae_assemblies_full.tsv`),
  downloaded (`data/references_raw/`), headers standardized to
  `{genome_id}__{contig}` (`data/references/`,
  mapping in `results/panel/contig_name_map.tsv`).
- `results/qc/assembly_stats.tsv` — `seqkit stats` for all 14, done.
- BUSCO (`sordariomycetes_odb10`) completed for all 14 genomes:
  consistently 97.9–98.2 % Complete — consistently high completeness, no
  outliers (details: `docs/decisions.md`).
- **Only 1 of 14 references has a supplied GFF3 annotation**
  (`GCA004346965_1`). The rest need the reannotation
  (BRAKER3/Liftoff) planned in Section 6.2 before gene-based
  orthogroup/PAV analyses (Section 7.1) become possible.

**Panel deliberately reduced to 5 genomes (deviation from the document, POC
scope decision, see `docs/decisions.md`):** The document assumes
all 14 genomes; for the POC, **one representative per
host type** is used instead, to reduce computational cost (repeat masking, annotation,
panel construction). Active `config/references.tsv` (5 rows):

| Host type | Representative | Rationale |
|---|---|---|
| Oryza | `7015` | canonical reference strain "70-15" |
| Triticum | `GCA036493215_1` (Br48) | explicitly declared T2T assembly |
| Wild grass | `LpKY97` | established reference strain from the literature |
| Eleusine | `GCA004346965_1` | only Complete-Genome option |
| Avena | `GCA059329645_1` | only Complete-Genome option |

The full 14-genome catalog remains archived in
`config/references_full_catalog_14genomes.tsv` (BUSCO results for
all 14 remain valid but are not carried further through the pipeline).

**Repeat masking (Section 6.3) — complete** for all 5
panel genomes (`envs/repeats.yaml`: RepeatModeler2 + RepeatMasker,
`workflow/rules/repeats.smk`; parallelized, ~2:48 h for the last 4
genomes). Genome-wide repeat fraction 10.1–16.5 %. **Four contigs identified as
strong mini-/accessory-chromosome candidates** (small + far
above-average repeat fraction): `LpKY97__CP050927.1` (3.0 Mb,
56.3 %), `LpKY97__CP050928.1` (0.9 Mb, 53.0 %),
`GCA059329645_1__CM181343.1` (1.3 Mb, 45.4 %),
`GCA059329645_1__CM181341.1` (1.2 Mb, 24.4 %) — details:
`docs/decisions.md`.

**Gene annotation (Section 6.2) — partial:** Liftoff transfers the
single real NCBI annotation (`GCA004346965_1`, 13,521 genes) onto all 5
panel genomes (94.6–99.98 % transferred successfully, loss correlates
sensibly with host distance). BRAKER3 de-novo annotation remains
blocked by the GeneMark license. **Confirmation of the mini-chromosome candidates:**
The four contigs previously flagged by repeat fraction also show
3–10× lower gene density than the genome average — repeat-rich AND
gene-poor, two independent lines of evidence. Details: `docs/decisions.md`.

**Phase II (Section 7) — partial:**
- **7.1 OrthoFinder:** done. 13,226 orthogroups, 92.7 % present in all 5
  genomes (strict_core). Panel prevalence thresholds (7.3) recalibrated
  for 5 instead of 14 references (`config/thresholds.yaml`).
- **7.2 Whole-genome alignment/SyRI:** done for 3 of 10 genome pairs
  (SyRI requires an identical chromosome count; LpKY97/GCA059329645_1 have
  more contigs than the other 3 genomes due to their presumed accessory
  chromosomes — nucmer/coverage still runs for all 10 pairs).
  Details: `docs/decisions.md`.
- **7.4 Starship/Captain candidates:** done (core result). `starfish
  annotate` found 68 HMM-validated YR/Captain genes across all 5 genomes.
  Synthesis with size/repeat context/cargo genes classifies 35 of these
  as `starship_like`. **Closing the loop:** The two LpKY97 mini-chromosome
  candidates previously found via repeat fraction + gene density
  (`CP050927.1`, `CP050928.1`) indeed carry Captain genes, several
  classified as `starship_like` — three independent lines of evidence
  converge on the same contigs. Details: `docs/decisions.md`.
- **7.3 Region types:** done. All 21,820 10-kb windows across the 5
  panel genomes classified (`results/panel/panel_regions.bed`):
  75.4% strict_core, 12.5% soft_core, 6.6% repeat_ambiguous, 1.6% shell,
  1.3% unclassified, 1.0% subtelomeric_dynamic, 1.0% accessory_chromosome,
  0.5% starship_like. `accessory_chromosome` windows lie exclusively
  on the previously found LpKY97/GCA059329645_1 mini-chromosome
  candidates — consistency check passed. Details: `docs/decisions.md`.

**Phase III (Section 8, panel construction) — done:**
`results/panel/Mo_multiref_panel_v1.fa` (+ `.fai`) and
`results/panel/panel_contig_manifest.tsv`: 3,677 deduplicated
panel regions (98% identity/90% mutual coverage, mmseqs2) from
4,213 pre-blocks — 1,783 strict_core, 1,610 soft_core, 245 shell, 30
starship_like, 5 accessory_chromosome, 4 private_accessory. **Important
finding:** The deduplication rate is lower than in the document's example, because
our panel deliberately contains maximally divergent host lineages rather
than close relatives (89% of strict_core clusters remain single-genome
entries, since >2% sequence divergence between host lineages is normal
even in core regions) — this makes the panel larger but more informative
for later cross-lineage mapping. Details: `docs/decisions.md`.

**Phase IV (Section 9, pilot-isolate selection) — done for 5/10:**
10 stratified pilot isolates selected (`config/samples.tsv`), 5 of
them with real downloaded/QC-checked raw data (B71, ZM12,
K23_123, E34, TF051MC7).

**Phase V (Section 10, long-read mapping) — done for the 5
available isolates:** All mapped against the panel (minimap2,
preset chosen per actual platform: map-pb for PacBio raw data, map-ont
for Nanopore). Mapping rates 71.6–99.5% (primary). **Most important
finding:** One tiny accessory region (a 35-kb contig from the
Avena reference genome) shows extremely high coverage (300–4,000×)
across all 5 isolates — despite different host lineages — a possible sign
of a cross-host multi-copy element (interpretation still
preliminary, MAPQ filtering still pending). Details: `docs/decisions.md`.

**Phase VI (Section 11, PAV analysis) — done, with a core result:**
`results/pav/pav_matrix.tsv` (3,677 panel regions × 5 isolates,
two evaluation tiers: MAPQ≥20-filtered vs. all alignments including
secondary, to detect multi-mapping artifacts). **399 regions are
present in all 5 isolates — including 2 `starship_like` regions and
the previously found Avena-specific accessory region.** Both
Starship candidates (30 kb each, high repeat content) are confirmed present
in all 5 test isolates from different host lineages
(Triticum, Eleusine, wild grass) — even after strict MAPQ filtering, so
not a pure multi-mapping artifact. **This positively demonstrates the
POC's core question (cross-host Starship detectability via long-read
mapping) with real data.** Details: `docs/decisions.md`.

**Phase VII (Section 12, SV calling) — done:** Sniffles2 per isolate
+ cohort merge. **193 structural variants** in the 5-isolate cohort
(110 deletions, 82 insertions, 1 inversion). Not yet merged with
the PAV matrix into combined evidence (Section 13).

**Phase IX (Section 14, rarefaction) — done, important finding:**
Exhaustive C(5,k) combination analysis (instead of 1000 random
permutations as in the document — stricter possible with only 5 genomes).
**The 5-genome panel is NOT saturated** — the gain from the 5th genome is
around ~19-20% for all candidate classes (well above the 2-5% threshold).
Confirmed by an unmapped-read analysis: B71 has 28% unmapped reads
(273 Mb, ~6× genome size) — consistent with the already observed
lower mapping rate. **Consequence:** A larger reference panel (e.g. from the
archived 14-genome catalog) would likely reveal substantially more
candidate regions. Details: `docs/decisions.md`.

**Phase VIII (Section 13, candidate-region manifest) — done:**
`results/panel/panel_candidate_regions.tsv` (39 candidate regions: 30
starship_like, 5 accessory_chromosome, 4 private_accessory, all with
stable IDs STAR_/ACC_) and `results/pav/candidate_region_calls.tsv`
(final per-isolate assignment, combining PAV + SV evidence). **All 30
starship_like candidates have an assigned Captain gene (100%
consistency between the Starfish finding and panel clustering).** Section
13.3 (gene functional annotation per PAV block) not implemented — a
panel-wide gene GFF3 is missing. Details: `docs/decisions.md`.

This means all core phases of the workflow document (I–IX) have been addressed.

## Critical finding on test-isolate selection (Section 9)

The document assumes "42 chromosome-level long-read isolates" as the
test-isolate pool. The actual NCBI catalog yields exactly 42
assemblies at chromosome level — but:

- **Only 15 of 42 are actually long-read/hybrid sequenced**
  (13 long-read + 2 hybrid); 25 are **short-read**-based (presumably
  reference-guided scaffolded, mostly a large batch of Brazilian
  wheat isolates), 2 are historical Sanger (70-15 duplicate
  GCA/GCF_000002495.2).
- **Host diversity in this pool is strongly skewed:** 31/42
  *Triticum*, only 5 *Oryza*, 1 each of *Lolium*/*Setaria*, **0 *Eleusine*** —
  the stratification required by the document in Section 9.1 ("2
  Eleusine-associated isolates") is NOT achievable from this pool.
- See `config/samples_candidate_pool.tsv` for the full list (now
  44 rows, the `platform` column shows the actual technology).

**Eleusine gap closed:** The separate raw-read SRA catalog
(`../data/ncbi_m_oryzae_sra_wgs_longread.tsv`, 193 runs) was searched
by BioSample host attributes — 2 real Eleusine isolates with
PacBio raw data were found (**K23/123**, Kenya, ≈74×; **E34**, Ethiopia,
≈189×), both added as new rows to `config/samples_candidate_pool.tsv`.
Important secondary finding: none of the original 42 chromosome-level
isolates has findable raw reads in this catalog —
K23/123 and E34 are currently the only pool entries with actually
loadable FASTQ (currently being downloaded and QC-checked,
`data/longreads/`). For the remaining host groups (Triticum, Oryza,
Lolium, Setaria), a targeted search for raw data must still be
carried out before the final 10-pilot-isolate selection — not yet
done, so `config/samples.tsv` still cannot be finally populated.

## Directory structure

Exactly per document Section 4:

```text
magnaporthe_multiref_pav/
├── config/            # config.yaml, references.tsv, samples_candidate_pool.tsv, thresholds.yaml
├── data/
│   ├── references/          # header-standardized working copies (data/references/{genome_id}.fa)
│   ├── references_raw/      # unmodified NCBI downloads (original FASTA, not committed)
│   ├── annotations/          # (empty, Phase I/II)
│   ├── annotations_raw/     # supplied GFF3, where available (currently: 1/14)
│   ├── longreads/            # (empty, Phase IV)
│   └── resources/
├── envs/               # core, annotation, wga, longreads, reporting
├── workflow/
│   ├── Snakefile
│   ├── rules/          # qc.smk (implemented), rest as documented stubs
│   └── scripts/        # rename_fasta_headers.py (implemented)
├── results/
└── logs/
```
