# Analyseplan: Starfish-, Starship- und Mini-Chromosomen-Vergleich

**Projekt:** Promotion – Barragan / *Magnaporthe oryzae*  
**Zielgruppe:** Einstieg ohne umfangreiche Bioinformatik-Vorerfahrung  
**Datum:** 31. August 2026

---

## 1. Überblick: Was wird gebaut?

Du entwickelst keine einzelne Analyse, sondern eine **reproduzierbare Pipeline**. Eine Pipeline ist eine fest dokumentierte Abfolge von Schritten: Jede Eingabedatei, jeder Parameter und jede erzeugte Ergebnisdatei ist definiert. Dadurch kann die Analyse später auf 1, 40 oder 417 Isolate angewendet werden, ohne Befehle manuell kopieren zu müssen.

Die Pipeline beantwortet drei verschiedene, aber miteinander verbundene Fragen:

| Analysebereich | Zentrale Frage | Typische Evidenz |
|---|---|---|
| Starfish / Starships | Gibt es YR-assoziierte große mobile Elemente? | YR-HMM-Treffer, Captain-Kandidaten, Genkontext, Cargo-Gene, Grenzen |
| Mini-Chromosomen | Welche Sequenzen sind wahrscheinlich akzessorische Mini-Chromosomen? | Größe, Repeat-Gehalt, Core-Gen-Armut, Coverage, Syntenie, Homologie |
| Vergleich/HGT | Welche Elemente sind zwischen Isolaten ungewöhnlich ähnlich oder geteilt? | PAV, Alignment-Coverage, Sequenzidentität, Core-vs-Element-Kontrast, phylogenetische Diskordanz |

Ein wichtiger Grundsatz lautet:

> Ein YR-Treffer ist nicht automatisch ein Starship. Ein kurzer Contig ist nicht automatisch ein Mini-Chromosom. Hohe Sequenzähnlichkeit ist nicht automatisch horizontaler Transfer.

Für belastbare Aussagen brauchst du mehrere unabhängige Evidenzlinien.

---

## 2. Zentrale Begriffe

### Assembly

Eine **Assembly** ist die rekonstruierte Genomsequenz eines Isolats. Sie besteht aus einer oder mehreren Sequenzen, die Contigs oder Chromosomen genannt werden.

### Contig

Ein **Contig** ist eine zusammenhängend assemb­lierte DNA-Sequenz. In einer sehr guten Long-read-Assembly kann ein Contig einem vollständigen Chromosom entsprechen. In einer fragmentierten Assembly kann ein Chromosom auf mehrere Contigs verteilt sein.

### GFF

Eine **GFF-Datei** (*General Feature Format*) beschreibt, wo Gene und andere Merkmale auf der DNA liegen. Sie enthält unter anderem:

- Sequenzname bzw. Chromosom/Contig
- Start- und Endposition eines Gens
- Strangrichtung
- Feature-Typ, z. B. `gene`, `mRNA`, `CDS`
- IDs und weitere Attribute

Die Contig-ID in der ersten GFF-Spalte muss exakt dem Header der zugehörigen FASTA entsprechen.

### YR und Captain-Gen

YR bedeutet **Tyrosin-Rekombinase**. Tyrosin-Rekombinasen können DNA-Rekombination vermitteln. Bei Starships ist ein YR-Gen häufig ein Captain-Gen, das mit der Mobilität des Elements verbunden ist. Nicht jedes YR-Gen ist zwingend Teil eines Starship.

### Starship

Ein **Starship** ist ein großes mobiles genetisches Element, das ein Captain-Gen und variabel zusätzliche Gene, oft als Cargo bezeichnet, tragen kann. Um ein Starship überzeugend zu annotieren, braucht man mehr als den Nachweis einer YR-Domäne: Genumgebung, mögliche Grenzen, Wiederholungen und eine plausible Elementstruktur sind wichtig.

### Mini-Chromosom

Ein **Mini-Chromosom** oder akzessorisches Chromosom ist eine chromosomale Einheit, die nicht zum konservierten Core-Genom gehört. Solche Sequenzen können sich stark zwischen Isolaten unterscheiden und Gene enthalten, die an Anpassung, Pathogenität oder horizontalem Transfer beteiligt sind.

---

## 3. Warum Snakemake?

**Snakemake** ist eine Workflow-Engine für Bioinformatik. Du definierst Regeln statt eine lange Liste manueller Befehle auszuführen.

Beispielidee:

```text
normalisierte FASTA
        ↓
Starfish YR-Annotation
        ↓
konsolidierte Annotation
        ↓
Starship-Kandidaten
        ↓
PAV- und Ähnlichkeitsmatrizen
```

Snakemake erkennt, welche Zwischenprodukte fehlen, und startet nur die nötigen Schritte neu. Es kann mehrere Isolate parallel verarbeiten und speichert pro Regel Logs.

### Vorteile

- Gleicher Ablauf für alle Isolate
- Kein manuelles Kopieren von Kommandos
- Wiederaufnahme nach Fehlern
- Parallele Ausführung
- Klare Abhängigkeiten
- Gute Grundlage für Supplement, Dissertation und Manuskript
- Exakte Dokumentation der verwendeten Dateien und Parameter

### Warum nicht ein großes Python-Skript?

Python ist sinnvoll für eigene Berechnungen und Tabellen. Ein einziges großes Skript wird jedoch schnell schwer wartbar: Es muss selbst prüfen, welche Dateien existieren, welche Schritte erneut laufen sollen und wie parallele Jobs koordiniert werden. Snakemake übernimmt genau diese Orchestrierung. Python wird als Hilfsmittel innerhalb einzelner Pipeline-Schritte eingesetzt.

---

## 4. Empfohlene Projektstruktur

Lege die Analyse als klar getrennte Verzeichnisse an:

```text
barragan/
├── config/
│   ├── samples.tsv
│   ├── paths.yaml
│   ├── parameters.yaml
│   └── references.yaml
├── workflow/
│   ├── Snakefile
│   ├── rules/
│   │   ├── qc.smk
│   │   ├── headers.smk
│   │   ├── starfish.smk
│   │   ├── annotation.smk
│   │   ├── minichromosomes.smk
│   │   ├── comparisons.smk
│   │   └── reporting.smk
│   ├── scripts/
│   │   ├── normalize_fasta_headers.py
│   │   ├── normalize_gff_seqids.py
│   │   ├── parse_starfish_results.py
│   │   ├── calculate_contig_metrics.py
│   │   ├── define_mchr_candidates.py
│   │   ├── build_pav_matrix.py
│   │   └── score_hgt_candidates.py
│   └── envs/
│       ├── starfish.yaml
│       ├── alignment.yaml
│       └── python.yaml
├── input/
│   ├── assemblies/
│   ├── annotations/
│   └── metadata/
├── results/
│   ├── qc/
│   ├── normalized/
│   ├── starfish/
│   ├── minichromosomes/
│   ├── comparisons/
│   └── figures/
├── logs/
├── docs/
│   ├── workflow.md
│   ├── decisions.md
│   └── provenance.md
└── README.md
```

### Was gehört wohin?

| Ort | Inhalt |
|---|---|
| `input/` | Rohdaten bzw. kontrollierte Eingaben; nicht still verändern |
| `config/` | Tabellen und Parameter, die den Lauf steuern |
| `workflow/` | Snakemake-Regeln, kleine Analyse-Skripte und Software-Umgebungen |
| `results/` | durch die Pipeline erzeugte finale und Zwischenresultate |
| `logs/` | vollständige Standardausgabe und Fehlermeldungen je Regel/Isolat |
| `docs/` | Entscheidungen, Methodenbeschreibung und Datenprovenienz |

---

## 5. Metadaten zuerst aufbauen

Bevor viele Analysen laufen, braucht jedes Isolat eine eindeutige, stabile Kennung. Erstelle `config/samples.tsv` mit mindestens:

```tsv
isolate_id	assembly_accession	assembly_fasta	gff	host	lineage	country	year	sequencing	include
ASM434696v1	GCA_004346965.1	/path/to/GCA_004346965.1_ASM434696v1_genomic.fna	/path/to/GCA_004346965.1_ASM434696v1_genomic.gff	unknown	unknown	unknown	unknown	NCBI	TRUE
```

### Bedeutung der wichtigsten Spalten

| Spalte | Zweck |
|---|---|
| `isolate_id` | Kurze, eindeutige ID für Tabellen, Abbildungen und Dateinamen |
| `assembly_accession` | NCBI-Herkunftskennung, z. B. `GCA_004346965.1` |
| `assembly_fasta` | Vollständiger Pfad zur Originalassembly |
| `gff` | Vollständiger Pfad zur Genannotation; leer, falls nicht verfügbar |
| `host`, `lineage`, `country`, `year` | Biologische Metadaten für spätere Interpretation |
| `sequencing` | z. B. Illumina, PacBio, Nanopore oder unbekannt |
| `include` | Steuerung, ob ein Isolat in einem bestimmten Lauf verwendet wird |

**Regel:** Eine ID darf nie nachträglich „still“ geändert werden. Wenn eine ID korrigiert wird, dokumentiere den Grund in `docs/provenance.md`.

---

## 6. Phase A: Eingaben prüfen und standardisieren

### A1. Assembly-Qualität erfassen

Zuerst wird pro Assembly eine Qualitätsübersicht erzeugt. Mindestens erfassen:

- Anzahl Contigs
- Gesamtgröße
- längster Contig
- N50
- Anzahl `N`-Basen
- GC-Gehalt
- verfügbare Annotation
- Sequenziertechnologie, wenn bekannt

**Warum?** Unterschiede in Contig-Zahl und Fragmentierung können spätere Unterschiede bei Mini-Chromosomen imitieren. Ein fehlendes mChr kann biologisch fehlen – oder schlicht nicht assembliert worden sein.

### A2. FASTA-Header eindeutig machen

Starfish erwartet bzw. profitiert von Headern im Schema:

```text
>genomeID_contigID
```

Beispiel:

```text
Original:      >CP034210.1
Normalisiert:  >GCA_004346965.1_CP034210.1
```

Dies ist für die Multi-Isolat-Analyse wichtig, weil Contig-IDs sonst zwischen Isolaten nicht eindeutig sind.

### A3. GFF synchron anpassen

Wenn die FASTA-Header geändert werden, muss die erste Spalte der GFF synchron geändert werden.

Beispiel:

```text
Originale GFF-Spalte 1:      CP034210.1
Normalisierte GFF-Spalte 1:  GCA_004346965.1_CP034210.1
```

**Qualitätskontrolle:** Nach der Anpassung müssen alle SeqIDs aus der GFF in der FASTA vorhanden sein. Andernfalls dürfen FASTA und GFF nicht gemeinsam verwendet werden.

### A4. Originaldaten unverändert lassen

Die NCBI-Dateien bleiben unverändert in einem Rohdatenordner. Die Pipeline schreibt normalisierte Kopien nach beispielsweise:

```text
results/normalized/{isolate_id}.fna
results/normalized/{isolate_id}.gff
```

Dadurch ist jederzeit nachvollziehbar, was verändert wurde.

---

## 7. Phase B: Starfish-Workflow

### B1. YR-Kandidaten finden

Der erste Starfish-Schritt sucht YR-Gene de novo. Dazu nutzt du:

```text
YRsuperfamRefs.faa
YRsuperfams.p1-512.hmm
```

Das Muster lautet:

```bash
starfish annotate \
  --assembly <assemblies.tsv> \
  --profile <YRsuperfams.p1-512.hmm> \
  --proteins <YRsuperfamRefs.faa> \
  --prefix <isolate>_YR \
  --idtag YR \
  --outdir <output_dir> \
  --tempdir <temp_dir> \
  --threads <n>
```

### Was geschieht intern?

1. MetaEuk sucht in der DNA nach möglichen proteincodierenden Genen, die den Referenz-YR-Proteinen ähneln.
2. `hmmsearch` prüft die vorhergesagten Proteine gegen das YR-HMM.
3. Treffer, die die HMM-Schwelle bestehen, werden als gefilterte Kandidaten ausgegeben.

### B2. Ergebnisstatus korrekt benennen

Verwende klare Statusnamen:

| Status | Definition |
|---|---|
| `YR_candidate` | HMM-validiertes YR-Gen aus `starfish annotate` |
| `contextual_Starship_candidate` | YR-Kandidat mit einem plausiblen genomischen Kontext |
| `bounded_Starship_candidate` | Kandidat mit zusätzlich plausiblen Grenzen/flankierenden Merkmalen |
| `high_confidence_Starship` | Vorab definierte Kriterien für Struktur und Evidenz erfüllt |

Diese Terminologie verhindert, dass ein vorläufiger Treffer überinterpretiert wird.

### B3. Existierende Genannotation einbeziehen

Falls NCBI-GFF-Dateien vorhanden und mit der normalisierten FASTA kompatibel sind, übergib sie über `-g`:

```text
genomeID<TAB>Pfad-zur-normalisierten-GFF
```

Warum? Die de-novo-YR-Suche findet primär Kandidaten für YR-Gene. Für Cargo-Gene und Nachbarschaftsanalysen brauchst du möglichst vollständige Genmodelle im gesamten Genom.

### B4. Nachfolgende Starfish-Schritte

Die lokale Starfish-Installation sollte mit `starfish --help` dokumentiert werden. Typische Schritte im Starfish-Workflow sind:

1. Annotation bzw. Konsolidierung der vorhandenen und neuen Genmodelle.
2. Bestimmung des Genkontexts rund um YR-Kandidaten.
3. Suche nach Einfügepositionen und/oder Grenzen.
4. Analyse flankierender Sequenzen und möglicher Wiederholungen.
5. Zusammenfassung der Kandidaten und möglicher Cargo-Gene.

**Wichtig:** Die exakten Subcommands, Parameter und Dateiformate werden vor Umsetzung aus der lokal installierten Version dokumentiert. Nicht alle Starfish-Versionen haben dieselben Namen oder Optionen.

### B5. Starfish-Output pro Isolat

Empfohlene finale Tabelle:

```text
results/starfish/starship_candidates.tsv
```

Mögliche Spalten:

```tsv
candidate_id	isolate_id	contig_id	start	end	strand	yr_gene_id	yr_hmm_evalue	yr_hmm_score	context_class	boundary_class	cargo_gene_count	cargo_annotations	confidence	comments
```

---

## 8. Phase C: Mini-Chromosomen finden

Mini-Chromosomen müssen unabhängig von Starfish untersucht werden. Starfish ist kein allgemeines Tool zur sicheren mChr-Erkennung.

### C1. Pro Contig Metriken berechnen

Erstelle eine Tabelle pro Contig, z. B.:

```text
results/minichromosomes/contig_metrics.tsv
```

Empfohlene Spalten:

```tsv
isolate_id	contig_id	length_bp	gc_fraction	gene_count	core_gene_count	repeat_fraction	TE_fraction	secreted_protein_count	effector_candidate_count	median_depth	depth_ratio_to_core	mchr_score	classification
```

Nicht alle Spalten sind von Beginn an verfügbar. Beginne mit Länge, GC-Gehalt und Genanzahl; erweitere dann schrittweise.

### C2. Evidenz für mChr-Kandidaten

Ein Contig wird nicht über eine einzige Schwelle als mChr bezeichnet. Verwende eine Kombination von Merkmalen:

- geringe Größe im Verhältnis zu den Core-Chromosomen
- niedrige Dichte konservierter Core-Gene
- hoher Repeat- oder Transposon-Anteil
- hohe Dichte von Effektoren oder sekretieren Proteinen, falls biologisch plausibel
- auffällige Coverage bzw. Copy Number, falls Rohreads vorhanden sind
- geringe Syntenie zu Core-Chromosomen
- starke Homologie zu bekannten akzessorischen Regionen oder mChrs

### C3. Arbeitsklassen

Nutze konservative Klassen statt einer voreiligen Ja/Nein-Entscheidung:

| Klasse | Bedeutung |
|---|---|
| `core_like` | offenbar Teil eines konservierten Core-Chromosoms |
| `accessory_candidate` | zeigt einzelne oder mehrere akzessorische Eigenschaften |
| `mChr_candidate` | mehrere Evidenzlinien sprechen für ein Mini-Chromosom |
| `uncertain` | Assemblyfragmentierung oder Evidenzlage erlaubt keine sichere Einordnung |

### C4. Umgang mit fragmentierten Assemblies

Ein mChr kann in einer schlechten Assembly auf mehrere Contigs verteilt sein. Deshalb darf die Analyse nicht ausschließlich fragen: „Welcher einzelne Contig ist ein mChr?“

Zusätzlich sollten homologe Sequenzblöcke betrachtet werden. Mehrere kleine Contigs desselben Isolats können gemeinsam Teile desselben mChr repräsentieren, wenn sie auf dasselbe Referenz-mChr oder denselben mChr-Cluster aus verschiedenen Isolaten alignieren.

---

## 9. Phase D: Vergleich zwischen Isolaten

### D1. Core-Genom als Hintergrund bestimmen

Bevor HGT interpretiert wird, braucht man einen Vergleichsmaßstab: Wie ähnlich sind die Isolate im konservierten Genom?

Mögliche Ansätze:

- Single-copy-Orthologe bestimmen und alignieren
- konservierte Core-Regionen identifizieren
- pairwise Distanzen aus Core-SNPs oder Core-Alignments berechnen

Output-Beispiel:

```text
results/comparisons/core_distance.tsv
```

### D2. mChr- und Starship-Sequenzen vergleichen

Zunächst kann ein schneller Vorfilter ähnliche Sequenzen gruppieren. Anschließend folgt ein präziser Vergleich per Sequenzalignment.

Empfohlene Prinzipien:

- beide Alignierungsrichtungen berücksichtigen
- Query- und Target-Coverage speichern
- pro Paar die prozentuale Identität speichern
- orientierte und kollineare Blöcke dokumentieren
- nicht nur den besten kurzen Treffer verwenden

Output-Beispiel:

```text
results/comparisons/pairwise_element_alignments.tsv
```

Mögliche Spalten:

```tsv
query_isolate	query_element	target_isolate	target_element	query_coverage	target_coverage	aligned_bp	identity	orientation	synteny_class
```

### D3. Presence/absence-Variation (PAV)

PAV bedeutet: Ist ein Element oder ein homologer Block in einem Isolat vorhanden, abwesend oder unklar?

Beispiel einer PAV-Matrix:

| Elementcluster | Isolat_A | Isolat_B | Isolat_C |
|---|---:|---:|---:|
| mChr_cluster_001 | present | present | absent |
| Starship_cluster_014 | present | uncertain | absent |

Bei Fragmentierung ist `uncertain` besser als ein falsches `absent`.

### D4. Core-vs-Element-Kontrast

Eine zentrale HGT-Logik lautet:

> Wenn zwei Isolate im Core-Genom relativ unterschiedlich sind, aber ein mChr oder Starship fast identisch und großflächig geteilt wird, ist das ein Hinweis auf einen jüngeren Austausch oder Transfer des Elements.

Beispielhafte Darstellung:

```text
Isolat A vs. B
Core-Genom-Identität:                  97.5 %
Identität eines mChr-Clusters:         99.95 %
reziproke mChr-Abdeckung:              92 %
```

Das Beispiel ist keine feste Entscheidungsregel. Die Schwellen müssen anhand deiner Assemblies, deiner Isolatpopulation und nach Möglichkeit anhand bekannter positiver und negativer Vergleichsfälle kalibriert werden.

### D5. Phylogenetische Diskordanz

Ein zusätzliches starkes Signal ist, wenn die Abstammungsbeziehung eines Elements nicht zur Core-Genom-Verwandtschaft passt.

Beispiel:

- Im Core-Genom gruppieren Isolat A und B nicht eng zusammen.
- Ihre mChr-Sequenzen gruppieren aber eng und haben hohe, breite Sequenzähnlichkeit.

Das kann zu einem Transfer passen. Alternativen wie konservierte Selektion, Kontamination, Fehlassemblierung oder unvollständige Stichprobe müssen aber geprüft werden.

---

## 10. Minimaler erster Pipeline-Meilenstein

Baue zuerst keinen vollständigen Großworkflow. Erreiche nacheinander diese überprüfbaren Meilensteine.

### Meilenstein 1: Ein Isolat vollständig vorbereiten

- `samples.tsv` anlegen
- Original-FASTA und gegebenenfalls Original-GFF lokalisieren
- FASTA-Header normalisieren
- GFF-SeqIDs synchron normalisieren
- prüfen, dass alle GFF-SeqIDs in der FASTA existieren

**Erfolgskriterium:** Die normalisierte FASTA und GFF passen zusammen und sind dokumentiert.

### Meilenstein 2: Starfish stabil auf einem Isolat ausführen

- YR-Annotation mit den normalisierten Contig-IDs wiederholen
- Log auf Header-Warnungen und Toolfehler prüfen
- gefilterte YR-GFF kontrollieren
- Ergebnis in eine tabellarische Kandidatenliste überführen

**Erfolgskriterium:** Ein reproduzierbarer Starfish-Output mit eindeutigen Sequenz- und Gen-IDs.

### Meilenstein 3: Drei Isolate als Pilot analysieren

Wähle, wenn möglich:

- ein Isolat mit bekannter oder erwarteter mChr-Evidenz
- ein nah verwandtes Isolat
- ein weiter entferntes Isolat

Vergleiche Contig-Metriken, Starfish-Kandidaten und Sequenzhomologien.

**Erfolgskriterium:** Die Pipeline läuft mehrfach ohne manuelle Änderung der Befehle; Unterschiede und Grenzen der Daten werden sichtbar.

### Meilenstein 4: Auf alle Isolate skalieren

Erst wenn der Pilot plausibel ist:

- alle Isolate in `samples.tsv` aktivieren
- verfügbare Threads und Speicher realistisch einstellen
- pro Schritt Logs und Checks erzeugen
- Ergebnisse in Gesamttabellen zusammenführen

---

## 11. Qualitätskontrollen und Fehlerquellen

### Fehlende oder fragmentierte mChrs

Wenn eine Assembly nur aus kurzen Contigs besteht, kann ein vorhandenes mChr fragmentiert oder gar nicht vollständig assembliert sein. „Nicht gefunden“ ist dann nicht zwingend „biologisch abwesend“.

### Unterschiedliche Sequenzierqualität

Long-read-Assemblies sind für die Erkennung vollständiger mChrs, Wiederholungen und Elementgrenzen meist besser geeignet als stark fragmentierte Short-read-Assemblies. Assemblierungsqualität muss daher in jeder Vergleichstabelle mitgeführt werden.

### Contamination und Fehlzuordnung

Ungewöhnlich ähnliche Elemente zwischen weit entfernten Isolaten können HGT zeigen, aber auch durch Kontamination, Probenverwechslung oder falsch zugeordnete Contigs entstehen. Prüfe bei starken Kandidaten Read-Mapping, Abdeckung, Taxonomie der Gene und die Assembly-Herkunft.

### Überinterpretation eines YR-Treffers

Ein YR-HMM-Treffer ist der Beginn, nicht das Ende der Analyse. Bewahre die Trennung zwischen YR-Kandidat, Starship-Kandidat und hochkonfidentem Element in allen Tabellen und Figuren.

### Schwellenwerte nach Ergebnissichtung ändern

Notiere Kriterien und Parameter frühzeitig in `docs/decisions.md`. Wenn später Änderungen notwendig werden, versioniere sie mit Begründung und wiederhole betroffene Analysen.

---

## 12. Dokumentation während der Arbeit

Erstelle und pflege mindestens drei Dokumente.

### `README.md`

Enthält:

- Ziel der Pipeline
- kurze Installationsanleitung
- Beispielaufruf
- Ordnerübersicht
- erwartete Hauptoutputs

### `docs/workflow.md`

Enthält:

- grafische oder textuelle Reihenfolge der Regeln
- Input und Output je Regel
- verwendete Tools
- erwartete Laufzeit-/Ressourcenordnung

### `docs/decisions.md`

Enthält methodische Entscheidungen, z. B.:

```text
2026-08-31
Entscheidung: FASTA-Header werden als <assembly_accession>_<original_contig_id> normalisiert.
Begründung: Eindeutige Contig-IDs für Multi-Isolat-Analysen und Kompatibilität mit Starfish.
Folge: Zugehörige GFF-SeqIDs werden synchron angepasst.
```

Dokumentiere außerdem:

- verwendete Starfish-Version
- Datenbankdateien und HMM-Schwellen
- MetaEuk-Parameter
- Kriterien für mChr-Klassen
- Kriterien für Starship-Konfidenz
- Kriterien für geteilte Elemente und HGT-Priorisierung

---

## 13. Konkrete nächste Arbeiten

Die sinnvollste Reihenfolge ist:

1. **Rohdaten inventarisieren.** Für jedes Isolat FASTA, GFF, mögliche Rohreads und Metadaten in `samples.tsv` erfassen.
2. **Normalisierung programmieren.** Ein kleines Python-Skript erstellen, das FASTA-Header sicher umbenennt und eine Mapping-Tabelle `old_id → new_id` schreibt.
3. **GFF-Synchronisierung programmieren.** Ein zweites Skript, das die Mapping-Tabelle verwendet und die SeqIDs in Spalte 1 ersetzt.
4. **Validierung einbauen.** Eine Regel, die prüft, dass GFF-SeqIDs und FASTA-SeqIDs übereinstimmen.
5. **Starfish-Test wiederholen.** `GCA_004346965.1` mit normalisierter FASTA testen, um die bisherigen Header-Warnungen zu beseitigen.
6. **Starfish-Folgeworkflow lokal dokumentieren.** `starfish --help` und die Hilfe aller verwendeten Subcommands speichern; die exakten Schritte erst danach als Snakemake-Regeln implementieren.
7. **Contig-Metriken erzeugen.** Länge, GC-Gehalt und Genanzahl für alle Contigs berechnen; dies ist der Startpunkt der mChr-Analyse.
8. **Pilotvergleich planen.** Drei repräsentative Isolate auswählen und ihre mChr-Kandidaten sowie YR/Starship-Kandidaten vergleichend alignieren.
9. **Erst dann skalieren.** Nach erfolgreichem Pilotlauf die gesamte Isolatkollektion rechnen.

---

## 14. Praktische Leitlinie

Arbeite iterativ und halte jeden Schritt klein:

```text
1 Isolat → 3 Isolate → gesamter Datensatz
1 Regel   → Testoutput → Dokumentation → nächste Regel
```

Wenn ein Ergebnis überraschend ist, prüfe zuerst die Eingaben, Header, Assemblyqualität und Mapping-Abdeckung. Ein sauber dokumentiertes negatives oder unsicheres Ergebnis ist wissenschaftlich wertvoller als eine schnelle, aber nicht überprüfbare Behauptung.
