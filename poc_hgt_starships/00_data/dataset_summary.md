# NCBI *Pyricularia oryzae* Assembly-Datensatz — Zusammenfassung

Erhoben am 2026-09-02 via NCBI Datasets REST API v2 (`genome/taxon/318829`,
txid deckt sowohl "Pyricularia oryzae" als auch das Synonym "Magnaporthe
oryzae" ab). **605 Assemblies gesamt** (Stand: heute; NCBI wächst laufend
weiter). Volldatensatz: `ncbi_pyricularia_oryzae_assemblies_full.tsv`.

## 1. Speicheraufwand

| | |
|---|---|
| Summe Genomgroesse (unkomprimiert, alle 605) | 24.244.937.851 bp ≈ **24,2 GB** |
| ... gzip-komprimiert (typ. ~4x fuer DNA) | ≈ **6,1 GB** |
| Mittlere Genomgroesse | 40,1 Mb (Median 39,9 Mb) |
| Spannweite | 35,1 – 49,0 Mb |

Reines FASTA-Sequenzvolumen; Annotations-/Read-Rohdateien kommen additiv
dazu (siehe Abschnitt 5).

## 2. Host-Diversitaet

29 normalisierte Host-Kategorien (nach Tippfehler-/Synonym-Bereinigung,
z. B. "Eleucine"→"Eleusine", "Urochloe"/"Urochla"→"Urochloa",
"rice"/"Rice"→"Oryza sativa"), gruppiert nach Pflanzengattung:

| Gattung | n Assemblies |
|---|---|
| *Oryza* (Reis) | 275 |
| *Triticum* (Weizen) | 59 |
| *Eleusine* (Fingerhirse/Goosegrass) | 41 |
| *Urochloa* (Signalgras) | 19 |
| *Lolium* (Weidelgras) | 16 |
| *Setaria* (Hirse) | 7 |
| *Melinis* | 6 |
| *Bromus* | 4 |
| *Eragrostis* | 3 |
| *Stenotaphrum* | 2 |
| *Avena*, *Hordeum*, *Festuca*, *Echinochloa*, *Panicum*, *Leersia* | je 1 |
| *Zingiber* (Ingwer — vermutlich Fehlzuordnung/Ausreisser, separat pruefen) | 1 |
| ohne Host-Attribut | 164 (27 %) |

→ Volle Detailtabelle inkl. Rohbezeichnungen: `host_diversity_summary.tsv`.
27 % der Assemblies haben KEIN Host-Attribut im BioSample-Datensatz —
das begrenzt, wie viele fuer eine host-lineage-stratifizierte Panelauswahl
direkt nutzbar sind, ohne die Publikation/das BioProject einzeln zu pruefen.

## 3. Contig-/Scaffold-Groesse nach Assembly-Level

| Assembly-Level | n | Contig-N50 (Median) | Scaffold-N50 (Median) |
|---|---|---|---|
| Complete Genome | 14 | 6,37 Mb | 6,37 Mb |
| Chromosome | 42 | 0,15 Mb | 6,19 Mb |
| Contig | 151 | 0,16 Mb | 0,16 Mb |
| Scaffold | 398 | 0,06 Mb | 0,12 Mb |

Nur **14 "Complete Genome"** (2,3 %) + **42 "Chromosome"**-Level (6,9 %) —
zusammen 56 Assemblies (9,3 %) mit chromosomen-nahem Zusammenbau. Die
grosse Mehrheit (66 %, n=398) liegt nur auf Scaffold-Level.

**Explizit Telomer-zu-Telomer (T2T) geflaggt: 0.** Kein Assembly traegt
im Namen/Kommentar einen T2T-Vermerk. **Gapless chromosomen-Level**
(Contig-Zahl == Chromosomenzahl, ein pragmatisches Naeherungskriterium
fuer "T2T-artig", aber ohne bestaetigte Telomer-Repeats an beiden Enden)
trifft auf **14 Assemblies** zu — alle 14 "Complete Genome"-Eintraege,
ueberwiegend long-read-basiert (11 reines Long-Read, 3 Hybrid).

## 4. Datentyp (Sequenziertechnologie)

| Kategorie | n | Anteil |
|---|---|---|
| Short-read (Illumina/BGISEQ/DNBSEQ) | 508 | 84,0 % |
| Long-read (PacBio/Nanopore) | 70 | 11,6 % |
| Hybrid (Long+Short) | 24 | 4,0 % |
| Sanger (historisch) | 2 | 0,3 % |
| Hybrid (Short+Sanger) | 1 | 0,2 % |

Kreuztabelle Assembly-Level × Technologie zeigt den erwarteten
Zusammenhang: fast alle Scaffold-Level-Assemblies (395/398) sind
Short-Read-basiert; Complete-Genome-Level ist ausschliesslich
Long-Read/Hybrid (14/14).

## 5. Zusaetzlich verfuegbare Rohdaten (bereits katalogisiert, aus frueherer Sitzung)

Neben den 605 Assemblies gibt es einen separaten, bereits vorhandenen
SRA-Rohdaten-Katalog (`data/ncbi_m_oryzae_sra_runs.tsv`, 3.902
Sequenzierlaeufe aller Bibliotheksstrategien):
- 1.754 WGS-Kurzread-Laeufe (Illumina) — `data/ncbi_m_oryzae_sra_wgs_illumina.tsv`
- 193 WGS-Langread-Laeufe (Nanopore/PacBio) — `data/ncbi_m_oryzae_sra_wgs_longread.tsv`
- 1.204 RNA-Seq, 435 ChIP-Seq, weitere Bibliothekstypen (siehe
  `docs/decisions.md`, Eintrag "Umfassender NCBI-Sequenzdatensatz")

Diese SRA-Laeufe sind grossteils NICHT dieselben Isolate wie die 605
Assemblies (Assemblies sind bereits zusammengebaute Genome, SRA-Laeufe
sind ueberwiegend Rohreads ohne eigene Assembly) — fuer Phase 1 relevant
als Kandidatenpool fuer zusaetzliche Isolate ohne eigenen Assembly-Eintrag.

## Dateien in diesem Verzeichnis

- `ncbi_pyricularia_oryzae_assemblies_full.tsv` — alle 605 Assemblies,
  27 Spalten (Accession, Organismus, Stamm, Assembly-Level, Sequenziertechnologie,
  Host, Isolationsquelle, Geo, Sammel­datum, Contig-/Scaffold-N50 etc.)
- `host_diversity_summary.tsv` — normalisierte Host-Haeufigkeitstabelle
