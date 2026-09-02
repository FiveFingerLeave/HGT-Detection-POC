# Phase 0 - Rohdaten

Rohgenome (Long-Read-Assemblies bevorzugt), NCBI-Downloads, Metadaten (Klonlinie, Wirtspflanze) fuer das POC-Panel (Ziel: 20-30 Isolate).

**Input:** NCBI-Suche (Assembly-Level Complete Genome/Chromosome, Technologie PacBio/Nanopore) - Ausgangspunkt sind die bereits vorhandenen Kataloge in `data/ncbi_m_oryzae_*.tsv` und die bereits heruntergeladenen Genome in `assemblies/`.
**Output:** kuratierte Isolatliste mit Metadaten (Klonlinie, Wirt, Assembly-Level, Technologie, Quelle).

**Status:** Vollständiger NCBI-Assembly-Katalog abgerufen (605 Genome,
2026-09-02) — siehe `dataset_summary.md` für Speicher-/Host-/
Assembly-Level-/Datentyp-Zusammenfassung, `ncbi_pyricularia_oryzae_assemblies_full.tsv`
für die Volldaten, `host_diversity_summary.tsv` für die normalisierte
Host-Häufigkeitstabelle. **Noch offen:** finale Panel-Kuration (Filterung
auf Complete-Genome-/Chromosome-Level + Klonlinien-/Host-Diversität, Ziel
20–30 Isolate).
