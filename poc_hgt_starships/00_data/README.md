# Phase 0 - Raw Data

Raw genomes (long-read assemblies preferred), NCBI downloads, metadata (clone lineage, host plant) for the POC panel (target: 20-30 isolates).

**Input:** NCBI search (assembly-level Complete Genome/Chromosome, technology PacBio/Nanopore) - starting point is the catalogs already available in `data/ncbi_m_oryzae_*.tsv` and the genomes already downloaded in `assemblies/`.
**Output:** curated isolate list with metadata (clone lineage, host, assembly level, technology, source).

**Status:** Full NCBI assembly catalog retrieved (605 genomes,
2026-09-02) — see `dataset_summary.md` for the storage/host/
assembly-level/data-type summary, `ncbi_pyricularia_oryzae_assemblies_full.tsv`
for the full data, `host_diversity_summary.tsv` for the normalized
host frequency table. **Still open:** final panel curation (filtering
to Complete-Genome/Chromosome level + clone-lineage/host diversity, target
20–30 isolates).
