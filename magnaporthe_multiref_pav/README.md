# Multireferenzpanel & Long-Read-PAV-Workflow

Implementierung von
[Dokumentation/multireferenzpanel_pav_workflow.md](../Dokumentation/multireferenzpanel_pav_workflow.md)
(19 Abschnitte: QC, Annotation, Pangenom-Klassifikation, Panel-Bau,
Long-Read-Mapping, fensterbasierte PAV, SV-Calling, Rarefaction).

## Status (2026-09-02)

**Phase I (Abschnitt 6.1, QC) — läuft/teilweise fertig:**
- 14 Complete-Genome-Referenzen aus dem NCBI-Gesamtkatalog
  (`../poc_hgt_starships/00_data/ncbi_pyricularia_oryzae_assemblies_full.tsv`)
  identifiziert, heruntergeladen (`data/references_raw/`), Header auf
  `{genome_id}__{contig}` standardisiert (`data/references/`,
  Zuordnung in `results/panel/contig_name_map.tsv`).
- `results/qc/assembly_stats.tsv` — `seqkit stats` für alle 14, fertig.
- BUSCO (`sordariomycetes_odb10`) läuft (`results/qc/{genome_id}_busco_summary.txt`
  je Genom) — rechenintensiv, kann pro Genom 10–30 Min dauern.
- **Nur 1 von 14 Referenzen hat eine mitgelieferte GFF3-Annotation**
  (`GCA004346965_1`). Die übrigen 13 brauchen die in Abschnitt 6.2
  vorgesehene Reannotation (BRAKER3/Liftoff), bevor genbasierte
  Orthogruppen-/PAV-Analysen (Abschnitt 7.1) möglich sind.

**Alle weiteren Phasen (Annotation, Repeat-Masking, Whole-genome-
Alignment/SyRI, Panel-Bau, Long-Read-Mapping, PAV, SV-Calling,
Rarefaction): noch nicht implementiert** — siehe die einzelnen
`workflow/rules/*.smk`-Dateien für den jeweiligen Statuskommentar und
die Voraussetzungen.

## Kritischer Befund zur Testisolat-Auswahl (Abschnitt 9)

Das Dokument geht von „42 chromosomenbasierten Long-Read-Isolaten" als
Testisolat-Pool aus. Der tatsächliche NCBI-Katalog liefert exakt 42
Assemblies auf Chromosome-Level — aber:

- **Nur 15 von 42 sind tatsächlich Long-Read-/Hybrid-sequenziert**
  (13 long-read + 2 hybrid); 25 sind **Short-Read**-basiert (vermutlich
  referenzgestützt gescaffoldet, v. a. eine große Charge brasilianischer
  Weizen-Isolate), 2 sind historisches Sanger (70-15-Duplikat
  GCA/GCF_000002495.2).
- **Host-Diversität in diesem Pool ist stark verzerrt:** 31/42
  *Triticum*, nur 5 *Oryza*, je 1 *Lolium*/*Setaria*, **0 *Eleusine*** —
  die vom Dokument in Abschnitt 9.1 geforderte Stratifizierung ("2
  Eleusine-assoziierte Isolate") ist aus diesem Pool NICHT erfüllbar.
- Siehe `config/samples_candidate_pool.tsv` für die volle Liste (42
  Zeilen, `platform`-Spalte zeigt die tatsächliche Technologie).

**Konsequenz:** Für echte Eleusine-Long-Read-Testisolate muss der
separate SRA-Rohdaten-Katalog
(`../data/ncbi_m_oryzae_sra_wgs_longread.tsv`, 193 Läufe) nach
BioSample-Host-Attributen durchsucht werden — noch nicht geschehen. Bis
dahin ist `config/samples.tsv` (die finalen 10 Pilotisolate) nicht
final befüllbar.

## Ordnerstruktur

Exakt nach Dokument-Abschnitt 4:

```text
magnaporthe_multiref_pav/
├── config/            # config.yaml, references.tsv, samples_candidate_pool.tsv, thresholds.yaml
├── data/
│   ├── references/          # header-standardisierte Arbeitskopien (data/references/{genome_id}.fa)
│   ├── references_raw/      # unveränderte NCBI-Downloads (Original-FASTA, nicht committen)
│   ├── annotations/          # (leer, Phase I/II)
│   ├── annotations_raw/     # mitgelieferte GFF3, wo vorhanden (aktuell: 1/14)
│   ├── longreads/            # (leer, Phase IV)
│   └── resources/
├── envs/               # core, annotation, wga, longreads, reporting
├── workflow/
│   ├── Snakefile
│   ├── rules/          # qc.smk (implementiert), Rest als dokumentierte Stubs
│   └── scripts/        # rename_fasta_headers.py (implementiert)
├── results/
└── logs/
```
