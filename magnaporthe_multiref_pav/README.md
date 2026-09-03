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

**Repeat-Masking (Abschnitt 6.3) — abgeschlossen** für alle 5
Panel-Genome (`envs/repeats.yaml`: RepeatModeler2 + RepeatMasker,
`workflow/rules/repeats.smk`; parallelisiert, ~2:48 h für die letzten 4
Genome). Genomweiter Repeat-Anteil 10,1–16,5 %. **Vier Contigs als
starke Mini-/Accessory-Chromosom-Kandidaten identifiziert** (klein +
weit überdurchschnittlicher Repeat-Anteil): `LpKY97__CP050927.1` (3,0 Mb,
56,3 %), `LpKY97__CP050928.1` (0,9 Mb, 53,0 %),
`GCA059329645_1__CM181343.1` (1,3 Mb, 45,4 %),
`GCA059329645_1__CM181341.1` (1,2 Mb, 24,4 %) — Details:
`docs/decisions.md`.

**Genannotation (Abschnitt 6.2) — teilweise:** Liftoff überträgt die
einzige echte NCBI-Annotation (`GCA004346965_1`, 13.521 Gene) auf alle 5
Panel-Genome (94,6–99,98 % erfolgreich übertragen, Verlust korreliert
sinnvoll mit Host-Distanz). BRAKER3-De-novo-Annotation bleibt an der
GeneMark-Lizenz blockiert. **Bestätigung der Mini-Chromosom-Kandidaten:**
Die vier zuvor per Repeat-Anteil geflaggten Contigs zeigen zusätzlich
3–10× niedrigere Gendichte als der Genomdurchschnitt — repeat-reich UND
genarm, zwei unabhängige Evidenzlinien. Details: `docs/decisions.md`.

**Phase II (Abschnitt 7) — teilweise:**
- **7.1 OrthoFinder:** fertig. 13.226 Orthogruppen, 92,7 % in allen 5
  Genomen (strict_core). Panel-Prävalenzschwellen (7.3) für 5 statt 14
  Referenzen neu kalibriert (`config/thresholds.yaml`).
- **7.2 Whole-genome-Alignment/SyRI:** fertig für 3 von 10 Genompaaren
  (SyRI verlangt gleiche Chromosomenzahl; LpKY97/GCA059329645_1 haben
  durch ihre vermuteten Accessory-Chromosomen mehr Contigs als die
  übrigen 3 Genome — nucmer/Coverage läuft trotzdem für alle 10 Paare).
  Details: `docs/decisions.md`.
- **7.4 Starship-/Captain-Kandidaten:** fertig (Kernergebnis). `starfish
  annotate` fand 68 HMM-validierte YR-/Captain-Gene über alle 5 Genome.
  Synthese mit Größe/Repeat-Kontext/Cargo-Genen klassifiziert 35 davon
  als `starship_like`. **Ringschluss:** Die beiden zuvor per Repeat-
  Anteil + Gendichte gefundenen LpKY97-Mini-Chromosom-Kandidaten
  (`CP050927.1`, `CP050928.1`) tragen tatsächlich Captain-Gene, mehrere
  als `starship_like` eingestuft — drei unabhängige Evidenzlinien
  konvergieren auf denselben Contigs. Details: `docs/decisions.md`.
- **7.3 Regionstypen:** fertig. Alle 21.820 10-kb-Fenster über die 5
  Panel-Genome klassifiziert (`results/panel/panel_regions.bed`):
  75,4% strict_core, 12,5% soft_core, 6,6% repeat_ambiguous, 1,6% shell,
  1,3% unclassified, 1,0% subtelomeric_dynamic, 1,0% accessory_chromosome,
  0,5% starship_like. `accessory_chromosome`-Fenster liegen ausschließlich
  auf den bereits gefundenen LpKY97-/GCA059329645_1-Mini-Chromosom-
  Kandidaten — Konsistenzprüfung bestanden. Details: `docs/decisions.md`.

**Phase III (Abschnitt 8, Panel-Bau) — fertig:**
`results/panel/Mo_multiref_panel_v1.fa` (+ `.fai`) und
`results/panel/panel_contig_manifest.tsv`: 3.677 deduplizierte
Panel-Regionen (98% Identität/90% gegenseitige Abdeckung, mmseqs2) aus
4.213 Vorab-Blöcken — 1.783 strict_core, 1.610 soft_core, 245 shell, 30
starship_like, 5 accessory_chromosome, 4 private_accessory. **Wichtiger
Befund:** Die Dedup-Rate ist niedriger als im Dokument-Beispiel, weil
unser Panel bewusst maximal divergente Host-Linien statt naher
Verwandter enthält (89% der strict_core-Cluster bleiben Einzelgenom-
Einträge, da >2% Sequenzdivergenz zwischen Wirtslinien selbst in
Core-Regionen normal ist) — macht das Panel größer, aber informativer
für späteres Cross-Lineage-Mapping. Details: `docs/decisions.md`.

**Phase IV (Abschnitt 9, Pilotisolat-Auswahl) — fertig für 5/10:**
10 stratifizierte Pilotisolate ausgewählt (`config/samples.tsv`), 5
davon mit real heruntergeladenen/QC-geprüften Rohdaten (B71, ZM12,
K23_123, E34, TF051MC7).

**Phase V (Abschnitt 10, Long-Read-Mapping) — fertig für die 5
verfügbaren Isolate:** Alle gegen das Panel gemappt (minimap2,
Preset je nach tatsächlicher Plattform: map-pb für PacBio-RAW, map-ont
für Nanopore). Mapping-Raten 71,6–99,5% (primär). **Wichtigster
Befund:** Eine winzige akzessorische Region (35-kb-Contig aus dem
Avena-Referenzgenom) zeigt bei allen 5 Isolaten — trotz unterschiedlicher
Wirtslinien — extrem hohe Coverage (300–4.000×), ein möglicher Hinweis
auf ein wirtsübergreifendes Multi-Kopie-Element (Interpretation noch
vorläufig, MAPQ-Filterung steht aus). Details: `docs/decisions.md`.

**Phase VI (Abschnitt 11, PAV-Analyse) — fertig, mit Kernergebnis:**
`results/pav/pav_matrix.tsv` (3.677 Panel-Regionen × 5 Isolate,
zwei Auswertungsebenen: MAPQ≥20-gefiltert vs. alle Alignments inkl.
secondary, um Multi-Mapping-Artefakte zu erkennen). **399 Regionen sind
bei allen 5 Isolaten präsent — darunter 2 `starship_like`-Regionen und
die zuvor gefundene Avena-spezifische Accessory-Region.** Beide
Starship-Kandidaten (je 30 kb, hoher Repeat-Anteil) sind bei allen 5
Testisolaten unterschiedlicher Wirtslinien (Triticum, Eleusine,
Wildgrass) bestätigt präsent — auch nach strenger MAPQ-Filterung, also
kein reines Multi-Mapping-Artefakt. **Damit ist die POC-Kernfrage
(wirtsübergreifende Starship-Nachweisbarkeit per Long-Read-Mapping) mit
echten Daten positiv demonstriert.** Details: `docs/decisions.md`.

**SV-Calling (Abschnitt 12), Rarefaction (Abschnitt 14): noch nicht
implementiert** — siehe die einzelnen `workflow/rules/*.smk`-Dateien für
den jeweiligen Statuskommentar.

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
