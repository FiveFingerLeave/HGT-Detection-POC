# Dokumentation: Vorbereitung des POC zur HGT- und Starship-Analyse in *Magnaporthe oryzae*

## Stand der Dokumentation

**Projekt:** Promotion – Populationsgenomik des horizontalen Gentransfers (HGT) bei *Magnaporthe oryzae*  
**Arbeitsphase:** Aufbau eines globalen, assembly-basierten Kandidatenpools für einen Proof of Concept (POC)  
**Arbeitsumgebung:** Windows Terminal mit Ubuntu über WSL; Miniforge/conda; isolierte Conda-Umgebung `starfish_env`  
**Dokumentationsstand:** Vor dem Download und der strukturellen Analyse ausgewählter Genome

---

## 1. Zielsetzung

Das Promotionsprojekt untersucht horizontalen Gentransfer (HGT) in *Magnaporthe oryzae* und *M. grisea*. Im Zentrum stehen zwei Klassen großer mobiler genomischer Elemente:

- akzessorische Mini-Chromosomen (mChr)
- Starships bzw. starship-ähnliche große mobile Elemente

Die übergeordnete Annahme lautet, dass HGT nicht rein zufällig geschieht, sondern aus wiederkehrenden Beziehungen zwischen Donorlinien, Empfängerlinien, mobilen Vektoren und lokalen ökologischen Bedingungen entsteht. Die spätere Analyse soll globale Donor–Empfänger–Vektor-Beziehungen identifizieren, deren lokale Wiederkehr in italienischen Feldpopulationen prüfen und HGT-Hotspots mit generalisierten linearen Mischmodellen (GLMM) modellieren.

Der aktuelle Arbeitsschritt ist ausdrücklich **kein vollständiger HGT-Nachweis und keine vollständige globale Analyse**. Es handelt sich um einen vorbereitenden Proof of Concept. Er beantwortet zunächst diese Machbarkeitsfrage:

> Gibt es genügend öffentlich verfügbare, ausreichend kontiguierliche und möglichst long-read-basierte *M. oryzae*-Assemblies, um Starships und mChr-Strukturen zuverlässig zu detektieren und daraus ein repräsentatives strukturelles Referenzpanel aufzubauen?

Die Beantwortung dieser Frage ist wichtig, weil Starship-Grenzen, Captain-Gene, flankierende Sequenzen, Insertionsstellen und mChr-Kontexte in stark fragmentierten Short-Read-Assemblies häufig nicht zuverlässig rekonstruiert werden können.

---

## 2. Methodischer Rahmen

### 2.1 Trennung zweier Analyseebenen

Der Workflow trennt bewusst eine strukturelle Assembly-Ebene von einer globalen Prävalenz-Ebene.

| Analyseebene | Zweck | Primäre Datengrundlage | Geeignete Methoden |
|---|---|---|---|
| Strukturelle Entdeckung | Vollständige bzw. weitgehend vollständige Starships, mChr-Kandidaten, Captains, Grenzen, flankierende Regionen und strukturelle Varianten erkennen | Long-Read- oder Hybrid-Assemblies mit hoher Contiguity | Starfish/Stargraph, Sequenzvergleich, Sytenie-/Alignments, mChr-Charakterisierung |
| Globale Prävalenz | Vorkommen, Presence/Absence und linienübergreifende Verteilung der zuvor definierten Kandidaten im breiten Isolatpanel messen | Short Reads bzw. bereits vorhandene Mapping-/Coverage-Daten | Mapping gegen Referenzen, Breadth-of-Coverage, PAV-Analyse, SNP-/Haplotypvergleich |

Die Trennung verhindert eine methodisch problematische Gleichsetzung von „in einer fragmentierten Assembly teilweise sichtbar“ und „strukturell als vollständiges mobiles Element belegt“.

### 2.2 POC-Logik

Der POC folgt einer Entscheidungslogik mit vier Modulen:

1. **Globale Frequenz-POC:** Identifikation möglicher identischer bzw. hochähnlicher mChr-/Starship-Elemente über verschiedene Linien hinweg.
2. **Lokale Rückfindungs-Approximation:** Abschätzung, wie viele der globalen Kandidaten in lokalen italienischen Feldpopulationen voraussichtlich wiedergefunden werden könnten.
3. **GLMM-Poweranalyse:** Simulation der erforderlichen Anzahl echter HGT-Ereignisse für ein stabiles Hotspot-Modell.
4. **Alternative Auswertungsstrategien:** Vereinfachung des Modells, Pooling von Vektorklassen oder explorative/deskriptive Analyse, falls die Ereigniszahl nicht genügt.

Der aktuelle Stand betrifft ausschließlich die technische und datenbezogene Vorbereitung von Modul 1.

---

## 3. Datengrundlage

### 3.1 Barragán-Supplementdaten

Als projektbezogene Ausgangsdatei liegt die Datei `Copy of Barragan2024_SupplementalTables.xlsx` vor. Relevante Sheets sind:

| Sheet | Inhalt | Nutzen im aktuellen Workflow |
|---|---|---|
| `TableS1` | Metadaten von neun italienischen klonalen Reisbrand-Isolaten | Lokaler Kontext: Isolat, Wirt, Ort, Jahr und Koordinaten |
| `TableS2` | Metadaten von 274 Isolaten mit genome-wide SNPs | Globaler Hintergrund und Strain-/Host-Metadaten |
| `TableS3` | Contig-Anzahl und Contig-Längen von neun italienischen Assemblies | Erste Einschätzung lokaler Assembly-Fragmentierung und mChr-Contigs |
| `TableS4` | Zusammenfassung der mChr-Contigs in den neun italienischen Isolaten | mChr-Kontext und bekannte Kandidaten |
| `TableS5` | Metadaten von 413 *M. oryzae*- und *M. grisea*-Isolaten | Zentrale Ausgangstabelle für spätere Verknüpfung: Strain ID, Wirt, Linie, Herkunft, BioProject und BioSample |
| `TableS6` | Genomweite und mChr-spezifische Breadth-of-Coverage-Werte über 413 Isolate | Beleg für die spätere globale coverage-/PAV-basierte Auswertung |
| `TableS10` | Long-Read-Sequenzierung und Assemblystatistiken für das *Eleusine*-Isolat Br62 | Dokumentierte Long-Read-Referenz mit Qualitätskennzahlen |

`TableS5` enthält insbesondere die Spalten `Strain ID`, `Host plant`, `M. oryzae Lineage`, `BioProject` und `Sample`. Die `Sample`-Spalte enthält BioSample-IDs wie `SAMN...` oder `SAMEA...`. Diese sind **keine Assembly-Accessions**, sondern Identifikatoren für biologische Proben und Metadaten.

### 3.2 NCBI-Assembly-Daten

Für den POC wurden keine FASTA-Dateien heruntergeladen. Stattdessen wurden zunächst nur NCBI-Assembly-Metadaten über das Taxon `Magnaporthe oryzae` abgefragt.

Ausgeführt wurde:

```bash
./datasets summary genome taxon "Magnaporthe oryzae" \
  --as-json-lines \
  > ncbi_m_oryzae_assemblies.jsonl
```

Ergebnis:

- **605 öffentliche NCBI-Assembly-Records**
- Datei: `ncbi_m_oryzae_assemblies.jsonl`
- Dateigröße: etwa **1,5 MB**
- Es wurden ausschließlich Metadaten geladen; keine Genome, Annotationen oder Rohreads.

Die anschließende TSV-Umwandlung erfasste folgende Felder:

```bash
./dataformat tsv genome \
  --inputfile ncbi_m_oryzae_assemblies.jsonl \
  --fields accession,organism-name,assminfo-name,assminfo-level,assminfo-assembly-method,assminfo-sequencing-tech,assminfo-bioproject,assminfo-biosample-accession,assminfo-release-date,assminfo-submitter,assmstats-total-sequence-len,assmstats-number-of-contigs,assmstats-contig-n50,assmstats-scaffold-n50,assmstats-genome-coverage \
  > ncbi_m_oryzae_assemblies.tsv
```

Damit entstand die Datei `ncbi_m_oryzae_assemblies.tsv`. Sie enthält eine Zeile pro NCBI-Assembly und bildet die Grundlage für die automatisierte Auswahl.

---

## 4. Verwendete Kennungen und Datenformate

| Präfix/Format | Bedeutung | Verwendung |
|---|---|---|
| `GCA_...` | GenBank-Assembly-Accession | Eindeutige Kennung einer konkret verfügbaren Genomassembly; korrekter Input für Assembly-Download |
| `GCF_...` | RefSeq-Assembly-Accession | NCBI-RefSeq-Version einer Assembly; nicht zusätzlich laden, wenn dieselbe biologische Assembly bereits als GCA gewählt ist |
| `SAMN_...` | NCBI BioSample | Verbindet Strain-/Probenmetadaten mit NCBI-Datensätzen; nicht direkt als Assembly-FASTA herunterladbar |
| `SAMEA_...` | ENA BioSample | Entsprechender BioSample-Typ im europäischen Archiv |
| `PRJNA_...` | NCBI BioProject | Projektcontainer mit vielen Proben, Reads und/oder Assemblies |
| `PRJEB_...` | ENA BioProject | Europäischer Projektcontainer |
| `SRR_...`, `ERR_...`, `DRR_...` | SRA-/ENA-/DDBJ-Run-Accessions | Rohreads; für den POC nicht primär herunterladen |
| `.fna` | FASTA-Genomdatei | Zentrale Eingabe für Starfish nach Auswahl geeigneter Assemblies |
| `.gff`/`.gff3` | Genomannotation | Optional für Gen-/Funktionskontext und Annotationen |
| `.gbff` | GenBank Flat File | Optionaler, umfassender Annotationscontainer |

---

## 5. Technische Umsetzung

### 5.1 Rechenumgebung

Die Arbeit erfolgt unter Windows mit Linux-Umgebung über WSL (Ubuntu). Die bioinformatische Umgebung liegt im Linux-Dateisystem und nicht im Windows-gemounteten Bereich unter `/mnt/c/...`.

Arbeitsverzeichnis:

```bash
~/promotion/barragan/data
```

Diese Ablage ist für bioinformatische Analysen vorteilhaft, da viele kleine Dateioperationen und spätere Mapping-/Assembly-Schritte im WSL-Dateisystem deutlich performanter als im gemounteten Windows-Dateisystem ausgeführt werden.

### 5.2 Softwareumgebung

Installiert und funktionsfähig:

- Miniforge als schlanke Conda-Distribution
- Conda-Environment: `starfish_env`
- Python 3.8 innerhalb dieser Umgebung
- Starfish als Tool zur Annotation großer mobiler Elemente bzw. Starship-ähnlicher Regionen
- NCBI Datasets CLI: `datasets`
- NCBI Metadata-Formatter: `dataformat`

Aktivierung der Arbeitsumgebung:

```bash
conda activate starfish_env
cd ~/promotion/barragan/data
```

Der erfolgreiche technische Test erfolgte über:

```bash
starfish --help
```

### 5.3 Assembly-Metadaten als TSV exportieren

Die verwendeten Metadatenfelder umfassen:

| Feld | Bedeutung im POC |
|---|---|
| `Assembly Accession` | Download-ID: GCA/GCF |
| `Organism Name` | Taxonomische Bezeichnung; kann *Pyricularia oryzae* statt *Magnaporthe oryzae* lauten |
| `Assembly Name` | Einreicher-/Assemblyname; kann einen Strain-Namen enthalten |
| `Assembly Level` | Grobe Qualitätsstufe: Contig, Scaffold, Chromosome oder Complete Genome |
| `Assembly Method` | Verwendeter Assembler, z. B. Canu, HGAP, NextDenovo, Celera oder SPAdes |
| `Assembly Sequencing Tech` | Sequenziertechnologie, z. B. Illumina, PacBio oder Oxford Nanopore |
| `Assembly BioProject Accession` | Projektverknüpfung |
| `Assembly BioSample Accession` | Zentraler Schlüssel für einen eindeutigen Merge mit `TableS5` |
| `Assembly Stats Total Sequence Length` | Plausibilitätscheck für Genomgröße |
| `Assembly Stats Number of Contigs` | Fragmentierungsmaß; niedriger ist besser |
| `Assembly Stats Contig N50` | Contiguity-Maß; höher ist besser |
| `Assembly Stats Scaffold N50` | Zusätzliches Contiguity-Maß |
| `Assembly Stats Genome Coverage` | Sequenziertiefe, sofern verfügbar |

---

## 6. Sequenziertechnologien im NCBI-Datensatz

Eine Frequenzanalyse der Technology-Spalte zeigte, dass der Gesamtdatensatz überwiegend aus Illumina-Assemblies besteht, aber einen substanziellen Long-Read-/Hybridanteil aufweist.

Beispiele der häufigsten gemeldeten Technologien:

| Technologie | Anzahl Records |
|---|---:|
| Illumina HiSeq | 134 |
| Illumina HiSeq 2000 | 116 |
| Illumina NovaSeq | 63 |
| Illumina HiSeq 2500 | 60 |
| Illumina MiSeq | 52 |
| Illumina GAIIx | 36 |
| Oxford Nanopore MinION | 27 |
| PacBio Sequel | 22 |
| Oxford Nanopore + Illumina | 12 |
| Oxford Nanopore (weitere Varianten) | mindestens 6 |
| PacBio (weitere Varianten) | mindestens 5 |
| PacBio-/Illumina-Hybridvarianten | mehrere weitere Records |

Die technischen Labels sind nicht vollständig standardisiert. Beispielsweise kommen `minION`, `MiniION`, `Oxford Nanopore MinION`, `PromethION`, `PacBio RSII`, `PacBio Sequel`, `PacBio; Illumina MiSeq` und weitere Varianten vor. Deshalb wird nicht nach einem exakten Begriff, sondern mit einem breiten, nicht case-sensitiven Muster gefiltert.

Verwendeter Filter:

```bash
grep -Ei 'pacbio|nanopore|minion|promethion|hybrid' \
  ncbi_m_oryzae_assemblies.tsv \
  > ncbi_m_oryzae_longread_candidates.tsv
```

Danach wurde die Header-Zeile ergänzt:

```bash
{
  head -n 1 ncbi_m_oryzae_assemblies.tsv
  cat ncbi_m_oryzae_longread_candidates.tsv
} > ncbi_m_oryzae_longread_candidates_with_header.tsv
```

Ergebnis:

- **95 Zeilen** in `ncbi_m_oryzae_longread_candidates_with_header.tsv`
- davon **1 Header + 94 Long-Read-/Hybrid-Kandidaten**

---

## 7. Qualitätskriterien für Starfish-Kandidaten

### 7.1 Begründung

Starships sind große mobile Elemente, häufig im Bereich von mehreren zehn bis mehreren hundert Kilobasen. Akzessorische Mini-Chromosomen können ebenfalls große, repeat-reiche und strukturell dynamische Sequenzräume umfassen. Stark fragmentierte Illumina-Assemblies mit tausenden Contigs und N50-Werten im Kilobasenbereich sind daher ungeeignet, um Elementgrenzen, strukturelle Variation und Integrationskontext belastbar zu bestimmen.

Ein beobachtetes Beispiel aus dem Metadatensatz sind Illumina-GAIIx-/Velvet-Assemblies mit rund 5.600–9.450 Contigs und Contig-N50-Werten von etwa 7,6–16 kb. Diese können für globale Mapping- oder PAV-Fragen weiterhin nützlich sein, bilden aber keine sinnvolle primäre Grundlage für die strukturelle Starship-Entdeckung.

### 7.2 Initialer, bewusst inklusiver Qualitätsfilter

Eine Long-Read-/Hybridassembly wurde in den ersten hochwertigen Kandidatenpool aufgenommen, wenn mindestens eines der folgenden Kriterien erfüllt war:

- `Assembly Level` entspricht `Chromosome` oder `Complete Genome`.
- Anzahl Contigs ist kleiner oder gleich 200.
- Contig-N50 ist mindestens 1.000.000 bp.

Der Filter lautet:

```bash
awk -F '\t' '
BEGIN { OFS="\t" }
NR == 1 { print; next }
$4 ~ /Chromosome|Complete Genome/ || ($12 != "" && $12 <= 200) || ($13 != "" && $13 >= 1000000) {
    print
}
' ncbi_m_oryzae_longread_candidates_with_header.tsv \
> ncbi_m_oryzae_longread_highquality.tsv
```

Die Kriterien sind absichtlich inklusiv, damit geeignete Long-Read-Assemblies nicht frühzeitig ausgeschlossen werden. Die endgültige Auswahl erfolgt danach durch Priorisierung und Kontrolle auf biologische Relevanz sowie mögliche Doppelungen.

### 7.3 Beobachtete hochwertige Beispiele

Die erste Sichtprüfung lieferte zahlreiche sehr hochwertige Kandidaten, etwa:

| Assembly | Technologie | Level | Contigs | Contig-N50 | Erste Priorität |
|---|---|---|---:|---:|---|
| `GCA_003015475.2` | PacBio Sequel | Contig | 13 | 5,50 Mb | High |
| `GCA_003015815.2` | PacBio Sequel | Contig | 13 | 5,45 Mb | High |
| `GCA_003015975.2` | PacBio Sequel | Contig | 10 | 5,64 Mb | High |
| `GCA_003016745.2` | PacBio Sequel | Contig | 16 | 6,16 Mb | High |
| `GCA_004346965.1` | PacBio RSII + Illumina | Complete Genome | 7 | 6,13 Mb | High |
| `GCA_004785725.2` | Oxford Nanopore | Chromosome | 10 | 6,47 Mb | High |
| `GCA_012272995.1` | Oxford Nanopore | Complete Genome | 9 | 6,16 Mb | High |
| `GCA_021442365.1` | PacBio + Illumina | Contig | 21 | 5,38 Mb | High |
| `GCA_021764705.1` | Oxford Nanopore + Illumina | Chromosome | 14 | 5,53 Mb | High |

Diese Beispiele zeigen, dass ein strukturelles Referenzpanel für den Starfish-POC technisch realisierbar ist.

### 7.4 Qualitätsklassen für die spätere Auswahl

| Kategorie | Kriterien | Verwendung |
|---|---|---|
| `high` | Long Read/Hybrid und Chromosome/Complete Genome, oder höchstens 30 Contigs, oder Contig-N50 mindestens 5 Mb | Primäre Starfish-Analyse; bevorzugte Referenzen |
| `medium` | Long Read/Hybrid und höchstens 200 Contigs oder Contig-N50 mindestens 1 Mb | Ergänzung des Panels; Einsatz nach Sichtprüfung |
| `exclude` | Keine klare Long-Read-/Hybridbasis oder unzureichende Contiguity | Nicht für primäre Strukturerkennung; ggf. später für Mapping/PAV |

---

## 8. Besondere Qualitätskontrollen

### 8.1 GCA/GCF-Doppelungen vermeiden

Eine GCA- und eine GCF-Accession können dieselbe zugrundeliegende biologische Assembly repräsentieren. Sie dürfen nicht als unabhängige Genome gezählt oder doppelt für Starfish heruntergeladen werden.

### 8.2 Wiederholte Einreichungen erkennen

Im ersten Kandidatenpool sind mehrere Assemblies mit identischen Kenngrößen sichtbar. Beispiele:

| Bevorzugte Assembly | Potenziell redundante Assembly | Indiz |
|---|---|---|
| `GCA_003015475.2` | `GCA_011799965.1` | gleiche Genomlänge, 13 Contigs, Contig-N50 5,50 Mb |
| `GCA_003015815.2` | `GCA_011799905.1` | gleiche Genomlänge, 13 Contigs, Contig-N50 5,45 Mb |
| `GCA_003015975.2` | `GCA_011799925.1` | gleiche Genomlänge, 10 Contigs, Contig-N50 5,64 Mb |
| `GCA_003016745.2` | `GCA_011799915.1` | gleiche Genomlänge, 16 Contigs, Contig-N50 6,16 Mb |

Diese Kandidaten müssen vor der finalen Auswahl über BioSample, Assembly-Name, Einreicherinformationen und ggf. NCBI-Assembly-Details kontrolliert werden. Solche Redundanzen dürfen nicht die Zahl biologisch unabhängiger Isolate erhöhen.

### 8.3 Taxonomische Synonyme berücksichtigen

Die NCBI-Metadaten können `Pyricularia oryzae` anzeigen, obwohl das Promotionsprojekt den Namen *Magnaporthe oryzae* verwendet. Das ist ein taxonomisches Synonym bzw. eine unterschiedliche Nomenklaturkonvention und muss beim Datenmerge nicht automatisch als Ausschlussgrund gewertet werden.

### 8.4 Long-Read-Label nicht blind vertrauen

Ein Long-Read-Technologieeintrag ist ein starkes Auswahlkriterium, aber kein vollständiger Qualitätsnachweis. Umgekehrt kann eine hochwertige Assembly trotz leerem Technology-Feld vorliegen. Deshalb erfolgt die Priorisierung immer kombiniert über:

- Sequenziertechnologie
- Assembly-Level
- Contig-Anzahl
- Contig-N50
- Scaffold-N50
- Genomlänge und Coverage als Plausibilitätswerte
- biologische Metadaten und Verknüpfung mit dem Barragán-Datensatz

---

## 9. Aktueller Datenbestand und Ordnerstruktur

Empfohlene bzw. verwendete Projektstruktur:

```text
~/promotion/barragan/
├── data/
│   ├── Copy of Barragan2024_SupplementalTables.xlsx
│   ├── ncbi_m_oryzae_assemblies.jsonl
│   ├── ncbi_m_oryzae_assemblies.tsv
│   ├── ncbi_m_oryzae_longread_candidates.tsv
│   ├── ncbi_m_oryzae_longread_candidates_with_header.tsv
│   └── ncbi_m_oryzae_longread_highquality.tsv
├── metadata/
│   ├── barragan_longread_master.xlsx
│   ├── assembly_accessions_selected.txt
│   └── download_manifest.tsv
├── assemblies/
│   ├── ncbi_dataset/
│   ├── genome_fasta/
│   ├── annotations/
│   └── sequence_reports/
├── starfish/
│   ├── input/
│   ├── output/
│   └── logs/
└── scripts/
    ├── build_master_metadata.py
    ├── select_assemblies.py
    └── run_starfish.sh
```

Bis zum aktuellen Stand sind nur Metadatendateien erzeugt worden. Es wurden keine FASTA-, GFF3-, Protein- oder FASTQ-Dateien gesammelt.

---

## 10. Nächste Arbeitsschritte

### Schritt 1: Master-Excel erzeugen

Der nächste unmittelbar geplante Schritt ist die Verknüpfung von:

- `TableS5` aus `Copy of Barragan2024_SupplementalTables.xlsx`, und
- `ncbi_m_oryzae_longread_highquality.tsv`.

Primärer Join-Key ist:

```text
TableS5: Sample  <->  NCBI: Assembly BioSample Accession
```

Die Master-Datei soll alle Barragán-Isolate behalten und pro Isolat den Match-Status dokumentieren:

| Status | Bedeutung |
|---|---|
| `matched` | Passende hochwertige Long-Read-/Hybridassembly über BioSample gefunden |
| `no_highquality_longread_match` | Kein Match im hochwertigen Long-Read-Pool; bedeutet nicht zwingend, dass keine Rohdaten oder Assembly existiert |
| `manual_check` | Unsicherer oder mehrdeutiger Abgleich; manuelle Prüfung nötig |

Zielname der Datei:

```text
barragan_longread_master.xlsx
```

### Schritt 2: Finales Starfish-Panel auswählen

Aus der Master-Datei wird ein erstes strukturelles Panel von etwa 20–40 biologisch unabhängigen Assemblies zusammengestellt. Auswahlkriterien:

- Repräsentation verschiedener Linien: Oryza, Triticum, *Eleusine*, Lolium und weitere Wildgraslinien.
- Höchste Qualität vor größter Anzahl.
- Keine GCA/GCF-Doppelungen und keine mehrfach eingereichten identischen Assemblies.
- Möglichst mehrere geographische und hostbezogene Kontexte.
- Einbezug von Isolaten mit Bezug zu Barragán-Studien, insbesondere bekannte mChrA-/mChr-Träger, sofern die Daten öffentlich verfügbar sind.
- Zusätzlich Referenzisolate wie 70-15, sofern eine hochwertige passende Assembly verfügbar ist.

### Schritt 3: Download der final ausgewählten Assemblies

Erst nach der Auswahl wird eine Datei `assembly_accessions_selected.txt` erzeugt, mit genau einer `GCA_...`- oder `GCF_...`-Accession pro Zeile.

Beispiel:

```text
GCA_003015475.2
GCA_004346965.1
GCA_004785725.2
```

Danach erfolgt der Download von FASTA, GFF3, Proteinen, GenBank-Flat-Files und Sequence-Reports:

```bash
./datasets download genome accession \
  --inputfile assembly_accessions_selected.txt \
  --include genome,gff3,protein,gbff,seq-report \
  --filename m_oryzae_starfish_panel.zip

mkdir -p ../assemblies/ncbi_dataset
unzip m_oryzae_starfish_panel.zip -d ../assemblies/ncbi_dataset
```

Die genomischen FASTA-Dateien können anschließend geprüft werden mit:

```bash
find ../assemblies/ncbi_dataset -type f -name "*.fna" | sort
```

### Schritt 4: FASTA-Dateien standardisieren und Qualitätskontrolle

Vor Starfish:

- Eindeutige Zuordnung jeder FASTA zu GCA/GCF, Strain-ID, BioSample und Linie.
- Prüfung der FASTA-Dateigröße, Contig-Anzahl und Header.
- Dokumentation von Ausschlüssen.
- Separate Ablage von Genom-FASTA, Annotation und Metadaten.

### Schritt 5: Starfish-Testlauf auf 1–3 Referenzassemblies

Zunächst keine Vollanalyse über das gesamte Panel. Empfohlen wird ein kontrollierter Testlauf auf:

- einer PacBio-Assembly mit sehr hoher Contiguity,
- einer ONT-/Hybridassembly,
- einer biologisch relevanten Referenz oder einem bekannten mChr-assoziierten Isolat.

Dabei wird geprüft:

- läuft Starfish technisch reproduzierbar;
- werden Captain-Kandidaten und starship-ähnliche Regionen ausgegeben;
- sind die Ergebnisse plausibel bezüglich Länge, Repeat-Kontext und Contiglage;
- welche Ausgabeformate für die anschließende Vergleichs- und Netzwerkphase entstehen.

### Schritt 6: Starship-Katalog und Sequenzvergleich

Nach erfolgreichem Test wird ein Kandidatenkatalog erzeugt, mindestens mit:

| Feld | Inhalt |
|---|---|
| `element_id` | Eindeutige Kandidatenkennung |
| `assembly_accession` | Herkunftsassembly |
| `strain_id` | Strain, soweit auflösbar |
| `lineage` | Wirt-/Linienzuordnung |
| `element_type` | Starship, starship-like, mChr-assoziiert, unklar |
| `contig` | Trägercontig |
| `start`, `end` | Koordinaten |
| `length_bp` | Länge |
| `captain_gene` | Captain-/Tyrosin-Rekombinase-Nachweis |
| `boundary_evidence` | z. B. TSD, terminale Sequenzen, flankierende Struktur |
| `quality_class` | High/Medium/Low |

Im nächsten Schritt werden die Elemente paarweise zwischen Linien verglichen. Für den POC ist ein pragmatisches Kriterium geeignet:

- mindestens 95% Nukleotididentität
- über mindestens 80% der Elementlänge
- Vorkommen in mindestens zwei unterschiedlichen hostassoziierten Linien

Diese Fälle sind noch nicht automatisch HGT-belegt. Sie sind **priorisierte HGT-Kandidaten**, die anschließend gegen die Kerngenom-Phylogenie, regionale Ähnlichkeit, PAV und ggf. D-Statistiken geprüft werden müssen.

### Schritt 7: Übergang zur globalen PAV-/Mapping-Ebene

Nur die strukturell gut definierten Kandidaten aus den hochwertigen Assemblies dienen anschließend als Referenzen für den großen globalen Short-Read-Datensatz. Dadurch wird die Frage beantwortbar:

> In welchen Linien, Hosts, Ländern und Isolaten sind dieselben oder hochähnliche mobilen Elemente präsent?

Dies ist der Punkt, an dem die bereits im Barragán-Datensatz vorhandenen Coverage-/mChr-Informationen (TableS6) und die globalen BioSample-/BioProject-Verknüpfungen besonders relevant werden.

---

## 11. Entscheidungskriterien für den POC

| POC-Frage | Positives Ergebnis | Konsequenz |
|---|---|---|
| Gibt es ausreichend hochwertige Assemblies? | Mindestens ca. 20 biologisch unabhängige Long-Read-/Hybridassemblies über mehrere Linien | Starfish-Panel erstellen und strukturelle Analyse starten |
| Gibt es starship-ähnliche/mChr-assoziierte Kandidaten? | Reproduzierbare Starfish-Kandidaten mit Captain-/Strukturmerkmalen | Kandidatenkatalog und Sequenzvergleich aufbauen |
| Gibt es linienübergreifend hochähnliche Elemente? | Elemente mit hoher Identität und hoher Abdeckung in mehreren Linien | Kandidaten in die HGT-Evidenzpipeline überführen |
| Sind die Ereignisse ausreichend häufig? | Mehrere unabhängige Donor-Empfänger-/Element-Beziehungen | Lokale Rückfindungs-Approximation und GLMM-Poweranalyse vorbereiten |
| Sind die Ereignisse zu selten? | Nur Einzelfälle bzw. keine wiederkehrenden Beziehungen | Modell vereinfachen; mChr/Starships poolen; explorative Ziele priorisieren; gezielte neue Proben prüfen |

---

## 12. Fazit zum aktuellen Stand

Der bisherige Workflow beantwortet die erste technische Machbarkeitsfrage positiv:

- Die NCBI-Abfrage liefert 605 öffentliche Assembly-Records für *Magnaporthe/Pyricularia oryzae*.
- Darunter wurden 94 Long-Read- oder Hybrid-Kandidaten anhand der Sequenziertechnologie identifiziert.
- Der Qualitätsfilter zeigt zahlreiche strukturell sehr geeignete Assemblies mit 7–56 Contigs und Contig-N50-Werten von etwa 3–6,5 Mb.
- Damit ist ein qualitativ belastbares, assembly-basiertes Starfish-Referenzpanel für einen POC realistisch.
- Der nächste kritische Arbeitsschritt ist nicht der vollständige Download, sondern die kontrollierte Verknüpfung dieser Assemblies mit den biologischen Metadaten aus Barragán TableS5 und die Entfernung biologischer bzw. technischer Doppelungen.

Die bisherige Arbeit ist damit als reproduzierbare Datengrundlagen- und Qualitätsphase einzuordnen. Ein HGT-Nachweis wurde noch nicht durchgeführt, und die beobachtete Assemblyverfügbarkeit darf noch nicht als Frequenz horizontaler Transfers interpretiert werden.
