# NCBI *Pyricularia oryzae* Assembly Dataset — Summary

Retrieved on 2026-09-02 via the NCBI Datasets REST API v2 (`genome/taxon/318829`,
this txid covers both "Pyricularia oryzae" and its synonym "Magnaporthe
oryzae"). **605 assemblies total** (as of today; NCBI keeps growing
continuously). Full dataset: `ncbi_pyricularia_oryzae_assemblies_full.tsv`.

## 1. Storage footprint

| | |
|---|---|
| Total genome size (uncompressed, all 605) | 24,244,937,851 bp ≈ **24.2 GB** |
| ... gzip-compressed (typ. ~4x for DNA) | ≈ **6.1 GB** |
| Mean genome size | 40.1 Mb (median 39.9 Mb) |
| Range | 35.1 – 49.0 Mb |

Pure FASTA sequence volume; annotation/raw-read files come on top
of this (see section 5).

## 2. Host diversity

29 normalized host categories (after cleaning typos/synonyms,
e.g. "Eleucine"→"Eleusine", "Urochloe"/"Urochla"→"Urochloa",
"rice"/"Rice"→"Oryza sativa"), grouped by plant genus:

| Genus | n assemblies |
|---|---|
| *Oryza* (rice) | 275 |
| *Triticum* (wheat) | 59 |
| *Eleusine* (finger millet/goosegrass) | 41 |
| *Urochloa* (signal grass) | 19 |
| *Lolium* (ryegrass) | 16 |
| *Setaria* (millet) | 7 |
| *Melinis* | 6 |
| *Bromus* | 4 |
| *Eragrostis* | 3 |
| *Stenotaphrum* | 2 |
| *Avena*, *Hordeum*, *Festuca*, *Echinochloa*, *Panicum*, *Leersia* | 1 each |
| *Zingiber* (ginger — likely misassignment/outlier, check separately) | 1 |
| no host attribute | 164 (27 %) |

→ Full detail table including raw labels: `host_diversity_summary.tsv`.
27 % of the assemblies have NO host attribute in the BioSample record —
this limits how many can be directly used for a host-lineage-stratified
panel selection without individually checking the publication/BioProject.

## 3. Contig/scaffold size by assembly level

| Assembly level | n | Contig N50 (median) | Scaffold N50 (median) |
|---|---|---|---|
| Complete Genome | 14 | 6.37 Mb | 6.37 Mb |
| Chromosome | 42 | 0.15 Mb | 6.19 Mb |
| Contig | 151 | 0.16 Mb | 0.16 Mb |
| Scaffold | 398 | 0.06 Mb | 0.12 Mb |

Only **14 "Complete Genome"** (2.3 %) + **42 "Chromosome"**-level (6.9 %) —
together 56 assemblies (9.3 %) with near-chromosome assembly. The
large majority (66 %, n=398) are only at scaffold level.

**Explicitly flagged telomere-to-telomere (T2T): 0.** No assembly carries
a T2T notation in its name/comment. **Gapless chromosome-level**
(contig count == chromosome count, a pragmatic approximation criterion
for "T2T-like", but without confirmed telomeric repeats at both ends)
applies to **14 assemblies** — all 14 "Complete Genome" entries,
predominantly long-read-based (11 pure long-read, 3 hybrid).

## 4. Data type (sequencing technology)

| Category | n | Share |
|---|---|---|
| Short-read (Illumina/BGISEQ/DNBSEQ) | 508 | 84.0 % |
| Long-read (PacBio/Nanopore) | 70 | 11.6 % |
| Hybrid (long+short) | 24 | 4.0 % |
| Sanger (historical) | 2 | 0.3 % |
| Hybrid (short+Sanger) | 1 | 0.2 % |

The assembly-level × technology cross-tabulation shows the expected
relationship: almost all scaffold-level assemblies (395/398) are
short-read-based; Complete-Genome level is exclusively
long-read/hybrid (14/14).

## 5. Additional raw data available (already cataloged, from an earlier session)

Besides the 605 assemblies, there is a separate, already existing
SRA raw-data catalog (`data/ncbi_m_oryzae_sra_runs.tsv`, 3,902
sequencing runs across all library strategies):
- 1,754 WGS short-read runs (Illumina) — `data/ncbi_m_oryzae_sra_wgs_illumina.tsv`
- 193 WGS long-read runs (Nanopore/PacBio) — `data/ncbi_m_oryzae_sra_wgs_longread.tsv`
- 1,204 RNA-Seq, 435 ChIP-Seq, and further library types (see
  `docs/decisions.md`, entry "Comprehensive NCBI sequence dataset")

These SRA runs are largely NOT the same isolates as the 605
assemblies (assemblies are already assembled genomes, SRA runs
are predominantly raw reads without their own assembly) — relevant for
Phase 1 as a candidate pool for additional isolates without their own assembly entry.

## Files in this directory

- `ncbi_pyricularia_oryzae_assemblies_full.tsv` — all 605 assemblies,
  27 columns (accession, organism, strain, assembly level, sequencing technology,
  host, isolation source, geo, collection date, contig/scaffold N50, etc.)
- `host_diversity_summary.tsv` — normalized host frequency table
