# Multireferenzpanel & Long-Read-PAV-Workflow

Implementierung von
[Dokumentation/multireferenzpanel_pav_workflow.md](../Dokumentation/multireferenzpanel_pav_workflow.md)
(19 Abschnitte: QC, Annotation, Pangenom-Klassifikation, Panel-Bau,
Long-Read-Mapping, fensterbasierte PAV, SV-Calling, Rarefaction).

## Status (2026-09-02)

**Phase I (Abschnitt 6.1, QC) — abgeschlossen für den vollen 14-Genom-Katalog:**
- 14 Complete-Genome-Referenzen aus dem NCBI-Gesamtkatalog
  (`../poc_hgt_starships/00_data/ncbi_pyricularia_oryzae_assemblies_full.tsv`)
  identifiziert, heruntergeladen (`data/references_raw/`), Header auf
  `{genome_id}__{contig}` standardisiert (`data/references/`,
  Zuordnung in `results/panel/contig_name_map.tsv`).
- `results/qc/assembly_stats.tsv` — `seqkit stats` für alle 14, fertig.
- BUSCO (`sordariomycetes_odb10`) für alle 14 Genome abgeschlossen:
  durchweg 97,9–98,2 % Complete — konsistent hohe Vollständigkeit, keine
  Ausreißer (Details: `docs/decisions.md`).
- **Nur 1 von 14 Referenzen hat eine mitgelieferte GFF3-Annotation**
  (`GCA004346965_1`). Die übrigen brauchen die in Abschnitt 6.2
  vorgesehene Reannotation (BRAKER3/Liftoff), bevor genbasierte
  Orthogruppen-/PAV-Analysen (Abschnitt 7.1) möglich sind.

**Panel bewusst auf 5 Genome reduziert (Abweichung vom Dokument, POC-
Scope-Entscheidung, siehe `docs/decisions.md`):** Das Dokument geht von
allen 14 Genomen aus; für den POC wird stattdessen **ein Repräsentant pro
Host-Typ** verwendet, um Rechenaufwand (Repeat-Masking, Annotation,
Panel-Bau) zu senken. Aktive `config/references.tsv` (5 Zeilen):

| Host-Typ | Repräsentant | Begründung |
|---|---|---|
| Oryza | `7015` | kanonischer Referenzstamm "70-15" |
| Triticum | `GCA036493215_1` (Br48) | explizit deklariertes T2T-Assembly |
| Wildgrass | `LpKY97` | etablierter Literatur-Referenzstamm |
| Eleusine | `GCA004346965_1` | einzige Complete-Genome-Option |
| Avena | `GCA059329645_1` | einzige Complete-Genome-Option |

Der volle 14-Genom-Katalog bleibt archiviert in
`config/references_full_catalog_14genomes.tsv` (BUSCO-Ergebnisse für
alle 14 bleiben gültig, werden aber nicht weiter durch die Pipeline
geführt).

**Repeat-Masking (Abschnitt 6.3) — läuft** für alle 5 Panel-Genome
parallel (`envs/repeats.yaml`: RepeatModeler2 + RepeatMasker,
`workflow/rules/repeats.smk`). Laufzeit pro Genom laut Pilotlauf
mehrere Stunden.

**Annotation, Whole-genome-Alignment/SyRI, Panel-Bau, Long-Read-Mapping,
PAV, SV-Calling, Rarefaction: noch nicht implementiert** — siehe die
einzelnen `workflow/rules/*.smk`-Dateien für den jeweiligen
Statuskommentar und die Voraussetzungen.

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
- Siehe `config/samples_candidate_pool.tsv` für die volle Liste (jetzt
  44 Zeilen, `platform`-Spalte zeigt die tatsächliche Technologie).

**Eleusine-Lücke geschlossen:** Der separate SRA-Rohdaten-Katalog
(`../data/ncbi_m_oryzae_sra_wgs_longread.tsv`, 193 Läufe) wurde nach
BioSample-Host-Attributen durchsucht — 2 echte Eleusine-Isolate mit
PacBio-Rohdaten gefunden (**K23/123**, Kenia, ≈74×; **E34**, Äthiopien,
≈189×), beide als neue Zeilen in `config/samples_candidate_pool.tsv`
ergänzt. Wichtiger Nebenbefund: keines der ursprünglichen 42
Chromosome-Level-Isolate hat auffindbare Rohreads in diesem Katalog —
K23/123 und E34 sind aktuell die einzigen Pool-Einträge mit real
ladbaren FASTQ (werden gerade heruntergeladen und QC-geprüft,
`data/longreads/`). Für die übrigen Host-Gruppen (Triticum, Oryza,
Lolium, Setaria) muss vor der finalen 10-Pilotisolate-Auswahl ebenfalls
gezielt nach Rohdaten gesucht werden — noch nicht geschehen, `config/
samples.tsv` ist deshalb weiterhin nicht final befüllbar.

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
