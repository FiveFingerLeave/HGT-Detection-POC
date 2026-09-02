# Reproduzierbarer Workflow: Multi-Referenzpanel und Long-Read-PAV-Analyse in *Magnaporthe oryzae*

## Ziel

Dieser Workflow baut aus **14 vollständigen Genomassemblies** ein versioniertes, biologisch annotiertes Multi-Referenzpanel. Es ermöglicht, Long Reads von Testisolaten gegen core-, accessory-, mini-chromosomale und Starship-like Genomregionen zu mappen, Presence/Absence Variation (PAV) zu bestimmen und jede Kandidatenregion einem oder mehreren Referenzisolaten zuzuordnen.

Anschließend werden zunächst **10 von 42 chromosomenbasierten Long-Read-Isolaten** aus unterschiedlichen Host-Typen als Pilotkohorte analysiert. Parallel wird eine Rarefaction-Analyse durchgeführt, um zu prüfen, ob das Referenzpanel die Akkumulation zusätzlicher Kandidatenregionen bereits sättigt.

---

## 1. Forschungsfragen

1. Welche Sequenzblöcke sind im 14-Genom-Set strict core, soft core, shell oder accessory?
2. Welche vollständigen oder partiellen accessory Chromosomen und Mini-Chromosomen kommen im Panel vor?
3. Welche Regionen sind Starship-like, also potenziell mobile, cargo-haltige Bereiche mit DUF3435-/Captain-Evidenz?
4. Welche dieser Regionen sind in den zehn Long-Read-Testisolaten vorhanden, abwesend oder aufgrund von Repeats/Mehrfachmapping nicht eindeutig bewertbar?
5. Welches bzw. welche Referenzisolate im Panel teilen eine in einem Testisolat nachgewiesene Kandidatenregion?
6. Führt das Hinzufügen weiterer Referenzgenome noch zu vielen neuen Kandidatenregionen oder nähert sich das Panel einer Sättigung?

---

## 2. Grundprinzip

Der Workflow trennt bewusst zwei Referenzprodukte:

1. **Vollständiger Referenzkatalog:** Alle 14 Genomassemblies bleiben vollständig, getrennt und unverändert erhalten. Er dient der Rekonstruktion von Syntenie, Herkunft, Homologie und Varianten.
2. **Analytisches Multi-Referenzpanel:** Ein nichtredundanter Satz von core-Repräsentanten plus alle validierten alternativen accessory Regionen, Mini-Chromosomen und Starship-like Regionen. Dieses Panel wird für Long-Read-Mapping und quantitative PAV-Analyse verwendet.

Nur alle 14 Genome direkt zu konkatenieren ist als exploratives Rohpanel möglich, aber nicht ideal für quantitative PAV-Calls: ähnliche core Regionen erzeugen Mehrfachmapping, verteilen Coverage auf mehrere Referenzen und können so falsche Abwesenheitscalls verursachen.

---

## 3. Inputs

| Datentyp | Umfang | Anforderung |
|---|---:|---|
| Vollständige Referenzassemblies | 14 Genome | chromosome-/contig-level; möglichst inklusive Mini-Chromosomen |
| Genannotation | 14 GFF3-Dateien und Proteome oder standardisierte Reannotation | Einheitliche IDs und konsistente Pipeline |
| Repeat-/TE-Annotation | 14 Genome | RepeatModeler2 + RepeatMasker oder EDTA |
| Long Reads | 42 Isolate verfügbar; 10 Pilotisolate | ONT oder PacBio HiFi; Metadaten erforderlich |
| Starship-Referenzen | DUF3435-HMM bzw. bekannte Captain-Proteine | Für Kandidatensuche und Annotation |
| Metadaten | Referenzen und Testisolate | Host, Hostgruppe, Herkunft, Jahr, Linie, Sequenztyp, Read-N50, Coverage |

---

## 4. Projektstruktur

```text
magnaporthe_multiref_pav/
├── config/
│   ├── config.yaml
│   ├── references.tsv
│   ├── samples.tsv
│   └── thresholds.yaml
├── data/
│   ├── references/
│   ├── annotations/
│   ├── longreads/
│   └── resources/
├── envs/
│   ├── core.yaml
│   ├── annotation.yaml
│   ├── wga.yaml
│   ├── longreads.yaml
│   └── reporting.yaml
├── workflow/
│   ├── Snakefile
│   ├── rules/
│   │   ├── qc.smk
│   │   ├── annotation.smk
│   │   ├── repeats.smk
│   │   ├── starships.smk
│   │   ├── wga.smk
│   │   ├── panel.smk
│   │   ├── mapping.smk
│   │   ├── pav.smk
│   │   ├── rarefaction.smk
│   │   └── report.smk
│   └── scripts/
│       ├── rename_fasta_headers.py
│       ├── classify_regions.py
│       ├── build_panel.py
│       ├── call_pav.py
│       ├── annotate_pav.py
│       └── rarefaction.py
├── results/
├── logs/
├── README.md
└── LICENSE
```

Snakemake ist als Workflow-Engine geeignet, weil es Regeln, reproduzierbare Software-Umgebungen, Konfigurationsdateien, parallele Ausführung, Logging und das Wiederaufnehmen unvollständiger Läufe unterstützt.

---

## 5. Metadaten und eindeutige IDs

### 5.1 Referenzmetadaten

Datei: `config/references.tsv`

```tsv
genome_id	host	host_group	lineage	country	year	assembly_fasta	annotation_gff	sequencing	complete_status
AG006	Oryza_sativa	oryza	Oryza_clonal_II	Italy	2018	data/references/AG006.fa	data/annotations/AG006.gff3	ONT_Illumina	complete
BR62	Eleusine_indica	eleusine	Eleusine	NA	NA	data/references/BR62.fa	data/annotations/BR62.gff3	ONT_Illumina	complete
...
```

### 5.2 Testisolatmetadaten

Datei: `config/samples.tsv`

```tsv
sample_id	host	host_group	lineage	country	year	platform	fastq	read_n50_bp	estimated_coverage	pilot
TEST01	Oryza_sativa	oryza	Oryza_clonal_I	Italy	2024	ONT	data/longreads/TEST01.fastq.gz	NA	NA	yes
TEST02	Triticum_aestivum	triticum	Triticum	NA	NA	ONT	data/longreads/TEST02.fastq.gz	NA	NA	yes
...
```

### 5.3 Globale FASTA-Header

Alle Contig-Namen müssen über sämtliche Referenzen eindeutig sein.

```text
>AG006__chr01
>AG006__chr02
>AG006__mChrA
>BR62__chr01
>BR62__contig07
```

Behalte die Original-FASTA-Dateien unverändert. Erzeuge neue, umbenannte Arbeitskopien und dokumentiere die Zuordnung in `results/panel/contig_name_map.tsv`.

---

## 6. Phase I: QC und Standardisierung der 14 Referenzen

### 6.1 Assembly-QC

```bash
seqkit stats data/references/*.fa > results/qc/assembly_stats.tsv

busco \
  -i data/references/AG006.fa \
  -l sordariomycetes_odb10 \
  -m genome \
  -o AG006_busco \
  -c 16
```

Dokumentiere pro Assembly:

- Gesamtgröße und erwartete Genomgröße
- Contig- und Chromosomenzahl
- Contig-N50 und größte Contigs
- BUSCO-Vollständigkeit, Fragmentierung und Duplikation
- Anteil unklassifizierter kleiner Contigs
- mitochondriale Contigs
- Repeat-Anteil
- telomerische Wiederholungen an Contig-Enden
- offensichtliche Kontamination
- bekannte bzw. vermutete Mini-Chromosomen

Mitochondrien, rDNA-Arrays und offensichtliche Kontaminationen werden nicht in das nukleäre PAV-Panel aufgenommen, aber separat archiviert.

### 6.2 Einheitliche Genannotation

Verwende möglichst eine standardisierte Genannotation für alle 14 Genome:

- BRAKER3 oder vergleichbare standardisierte ab-initio-/evidence-basierte Annotation
- Liftoff zur Übertragung hochwertiger Genmodelle auf nahe verwandte Genome
- InterProScan oder eggNOG-mapper zur Funktionsannotation
- SignalP und DeepTMHMM optional für sekretierte Proteine/Effektorkandidaten
- OrthoFinder für Orthogruppen

Wichtig: Lift-over allein kann accessory Gene unterschätzen. Ergänze daher durch de-novo-Annotation, besonders in nichtsyntenischen, subtelomerischen, mini-chromosomalen und repeat-reichen Bereichen.

### 6.3 Repeat- und TE-Annotation

Führe RepeatModeler2 + RepeatMasker oder EDTA konsistent auf allen Referenzen aus. Berechne den Repeat-Anteil pro Contig und pro Analysefenster.

Repeatreiche Bereiche werden nicht entfernt, aber als potenziell mehrdeutig markiert. Sie dürfen nicht allein auf Grundlage geringer eindeutiger Coverage als abwesend interpretiert werden.

---

## 7. Phase II: Pangenomische Klassifikation der Regionen

### 7.1 Orthogruppen

```bash
orthofinder \
  -f data/annotations/proteomes/ \
  -S diamond \
  -M msa \
  -T iqtree \
  -t 32 \
  -a 16
```

Nutze die Orthogruppenmatrix für eine genbasierte PAV-Klassifikation. Ein fehlendes Genmodell ist jedoch keine endgültige Evidenz für echten Verlust. Kandidaten werden daher zusätzlich DNA-basiert gegen die jeweiligen Assemblies validiert.

### 7.2 Whole-genome Alignments

Verwende mehrere Ankerreferenzen, mindestens je eine aus Oryza-, Triticum-, Eleusine- und Wildgras-assoziierten Linien. Ein einzelner Oryza-Anker würde die Diversität der übrigen Linien systematisch als fehlend oder nichtsyntenisch darstellen.

```bash
nucmer \
  --maxmatch \
  -l 100 \
  -c 500 \
  -p results/wga/AG006_vs_BR62 \
  data/references/AG006.fa \
  data/references/BR62.fa

delta-filter \
  -m \
  -i 85 \
  -l 500 \
  results/wga/AG006_vs_BR62.delta \
  > results/wga/AG006_vs_BR62.filtered.delta

show-coords \
  -THrd \
  results/wga/AG006_vs_BR62.filtered.delta \
  > results/wga/AG006_vs_BR62.coords.tsv
```

MUMmer4/nucmer ist für schnelle Whole-genome-Alignments geeignet; SyRI kann auf Whole-genome-Alignments syntenische Bereiche, Umordnungen und lokale Varianten unterscheiden. [web:37][web:44]

```bash
syri \
  -c results/wga/AG006_vs_BR62.coords.tsv \
  -r data/references/AG006.fa \
  -q data/references/BR62.fa \
  --prefix results/syri/AG006_vs_BR62_
```

### 7.3 Regionstypen

Jedes Intervall im Panel erhält eine Klasse in `results/panel/panel_regions.bed`.

| Regionstyp | Definition |
|---|---|
| `strict_core` | Homolog und weitgehend syntenisch in 14/14 Referenzgenomen |
| `soft_core` | Homolog in mindestens 13/14 Referenzgenomen |
| `shell` | In 3–12/14 Referenzgenomen vorhanden |
| `private_accessory` | Nur in einem oder zwei Referenzen vorhanden |
| `accessory_chromosome` | Ganzer Contig/Chromosom mit geringer Prävalenz, fehlender stabiler Kernsyntenie oder hoher struktureller Variabilität |
| `mini_chromosome` | Kleines eigenständiges Chromosom, idealerweise telomerbegrenzbar und klar vom Core getrennt |
| `subtelomeric_dynamic` | Region nahe Chromosomenende mit hoher Repeat- und PAV-Rate |
| `starship_like` | Mobiler Kandidatenblock mit DUF3435-/Captain-Evidenz, cargo-Genraum und variabler Präsenz |
| `repeat_ambiguous` | Bereich, dessen Präsenz wegen hoher Repeats oder Multi-Mapping nicht binär belastbar ist |
| `unclassified` | Noch nicht ausreichend klassifiziert |

Empfohlene Prävalenzschwellen bei 14 Referenzen:

```text
strict core: 14/14
soft core:   13/14
shell:       3–12/14
accessory:   1–2/14
```

### 7.4 Starship-like Kandidaten

Eine Region wird nur als `starship_like` geführt, wenn mehrere Evidenzlinien zusammenkommen:

- DUF3435-/Captain-Kandidat innerhalb oder nahe einer mutmaßlichen Grenze
- ausreichend großer, zusammenhängender Bereich, zum Beispiel mindestens 20 kb
- cargo-Gene bzw. funktionell heterogene Geninhalte
- fehlende stabile Core-Syntenie oder stark variable Präsenz
- optional terminale Wiederholungen oder konservierte Boundaries
- repeat-/TE- und/oder insertionreiche Umgebung

Ein DUF3435-Treffer allein reicht nicht. Bis zur strukturellen Bestätigung lautet die konservative Bezeichnung stets **Starship-like**.

---

## 8. Phase III: Aufbau des Multi-Referenzpanels

### 8.1 Zwei Panel-Ebenen

```text
A. Vollständiger Katalog:
   14 vollständige, getrennte Assemblies.

B. Analytisches Panel:
   - ein Repräsentant je strict-core-Homologieblock
   - alternative soft-core-/shell-Blöcke, falls nicht nahezu identisch
   - alle validierten accessory Regionen
   - alle Mini-Chromosomen
   - alle Starship-like Regionen
```

### 8.2 Deduplikationsregeln

- Strict-core-Homologieblöcke: eine Repräsentantensequenz wählen.
- Soft-core-/shell-Blöcke: alternative Varianten behalten, wenn sie biologisch oder strukturell verschieden sind.
- Accessory, Mini-Chromosomen und Starship-like Regionen: sämtliche validierten, nichtredundanten Varianten behalten.
- Hochidentische Sequenzen vor Aufnahme clustern, beispielsweise bei mindestens 98 % Identität und mindestens 90 % gegenseitiger Abdeckung.
- Bei homologen Regionen alle tragenden Referenzisolate als Mitglieder eines Regionclusters speichern.

### 8.3 FASTA-IDs

```text
>PANEL000104|type=mini_chromosome|cluster=ACC_018|rep=AG006|source=AG006__mChrA
>PANEL000105|type=starship_like|cluster=STAR_004|rep=BR62|source=BR62__chr05:2210000-2460000
>PANEL000106|type=strict_core|cluster=CORE_0321|rep=OryzaRef1|source=OryzaRef1__chr03:100000-450000
```

### 8.4 Panel-Manifest

Datei: `results/panel/panel_contig_manifest.tsv`

```tsv
panel_id	region_cluster	region_type	representative_reference	source_interval	length_bp	n_reference_genomes	reference_isolates	host_groups	repeat_fraction	starship_evidence
PANEL000104	ACC_018	mini_chromosome	AG006	AG006__mChrA:1-1200000	1200000	2	AG006;BR62	oryza;eleusine	0.58	no
PANEL000105	STAR_004	starship_like	BR62	BR62__chr05:2210000-2460000	250000	3	BR62;X12;X14	eleusine;wildgrass	0.43	DUF3435+boundary+PAV
PANEL000106	CORE_0321	strict_core	OryzaRef1	OryzaRef1__chr03:100000-450000	350000	14	ALL	all	0.08	no
```

Diese Datei ist der Schlüssel zur späteren biologischen Interpretation: Sie verbindet Panelkoordinaten mit Regionstyp, ursprünglicher Referenzkoordinate, Prävalenz, Repeat-Anteil und allen Referenzisolaten, die einen homologen Block tragen.

### 8.5 Finale Panel-FASTA

```bash
cat results/panel/panel_core.fa \
    results/panel/panel_accessory.fa \
  > results/panel/Mo_multiref_panel_v1.fa

samtools faidx results/panel/Mo_multiref_panel_v1.fa
```

Versioniere jede Änderung am Panel, zum Beispiel `v1.0`, `v1.1`, `v2.0`. Ändere alte Panelversionen nicht nachträglich.

---

## 9. Phase IV: Auswahl und QC der zehn Pilotisolaten

### 9.1 Stratifizierte Auswahl

Die zehn Isolate sollen die Diversität der 42 verfügbaren Long-Read-Isolate möglichst gut abbilden:

- 2–3 Oryza-assoziierte Isolate aus unterschiedlichen Klonen/Populationen
- 2 Triticum-assoziierte Isolate
- 2 Eleusine-assoziierte Isolate
- 2–3 Wildgras- bzw. weitere Host-Isolate
- mindestens ein Isolat mit erwarteter accessory DNA oder Mini-Chromosomen
- mindestens ein Benchmark-Isolat mit vorhandener vollständiger Assembly, falls möglich

### 9.2 Long-Read-QC

```bash
NanoPlot \
  --fastq data/longreads/TEST01.fastq.gz \
  --outdir results/qc/TEST01_nanoplot \
  --threads 8

seqkit stats data/longreads/TEST01.fastq.gz \
  > results/qc/TEST01_read_stats.tsv
```

Erfasse pro Isolat:

- Gesamtbasen und geschätzte Coverage
- Read-N50 und mittlere Readlänge
- Qualitätsverteilung
- Anteil sehr kurzer Reads
- erwartete nukleäre versus mitochondriale Coverage

Für belastbare Abwesenheitscalls sollte die gesamte nukleäre Coverage ausreichend sein; ein sinnvoller ONT-Startwert für den Pilotversuch ist etwa 20×, muss aber anhand deiner Daten kalibriert werden.

---

## 10. Phase V: Mapping der Long Reads

### 10.1 ONT-Mapping

```bash
minimap2 \
  -ax map-ont \
  --secondary=yes \
  -t 32 \
  results/panel/Mo_multiref_panel_v1.fa \
  data/longreads/TEST01.fastq.gz \
| samtools sort -@ 8 \
  -o results/mapping/TEST01.panel.bam

samtools index results/mapping/TEST01.panel.bam

samtools flagstat results/mapping/TEST01.panel.bam \
  > results/mapping/TEST01.flagstat.txt

samtools coverage results/mapping/TEST01.panel.bam \
  > results/mapping/TEST01.coverage_by_contig.tsv
```

Für PacBio HiFi wird `-ax map-hifi` verwendet.

### 10.2 Zwei Auswertungsebenen

1. **Alle Alignments einschließlich secondary Alignments:** Nachweis geteilter Homologie und Bewertung von Referenzalternativen.
2. **Primäre, hochwertige, eindeutige Alignments:** Quantitative PAV-Bewertung, zum Beispiel MAPQ ≥ 20 oder 30.

Eine Region mit hoher Homologie zu mehreren Referenzen darf nicht allein über den primären Alignmenttreffer einer einzelnen Referenz zugeschrieben werden.

---

## 11. Phase VI: Fensterbasierte PAV-Analyse

### 11.1 Fenster erzeugen

```bash
cut -f1,2 results/panel/Mo_multiref_panel_v1.fa.fai \
  > results/panel/panel.genome

bedtools makewindows \
  -g results/panel/panel.genome \
  -w 10000 \
  > results/panel/panel_10kb_windows.bed
```

Verwende 10-kb-Fenster als Standard für genomeweite PAV-Karten. Ergänze bei Starship-Grenzen, kleinen Genclustern oder Breakpoint-Validierung eine höhere Auflösung mit 1–2-kb-Fenstern.

### 11.2 Coverage berechnen

```bash
mosdepth \
  --by results/panel/panel_10kb_windows.bed \
  --threads 16 \
  results/mapping/TEST01 \
  results/mapping/TEST01.panel.bam
```

Pro Fenster werden berechnet:

- mittlere und mediane Tiefe
- Breadth bei mindestens 1×, 3× und 5×
- Anteil hochwertiger Alignments
- Anteil primärer versus sekundärer Alignments
- Repeat-Anteil
- Zahl eindeutiger Fenster bzw. diagnostischer k-mer
- Zahl der Reads, die Fenstergrenzen überspannen

### 11.3 PAV-Regeln

| Regionstyp | Present | Absent | Zusatzregel |
|---|---|---|---|
| Strict-/soft-core | Breadth ≥ 0,80 bei mindestens 5× | Breadth < 0,10 | Nur bei ausreichender globaler Coverage |
| Accessory-Block | Breadth ≥ 0,70 in mindestens 80 % der einzigartigen Fenster | Breadth < 0,10 in mindestens 80 % der einzigartigen Fenster | Multi-Mapping explizit auswerten |
| Mini-Chromosom | Stützung in mehreren entlang der Länge verteilten einzigartigen Fenstern | Keine stabile Coverage über Länge und eindeutige Absenz der Marker | Nicht aus einzelnen TE-Fenstern ableiten |
| Starship-like | Innenbereich plus mindestens eine Boundary unterstützt | Innenbereich und Boundaries nicht unterstützt | Boundary-Support erhöht Konfidenz |
| Repeat-ambiguous | Kein binärer Call | Kein binärer Call | Status `ambiguous` |

Benachbarte Fenster mit gleichem Status werden zu PAV-Blöcken zusammengeführt. Eine projektweite Mindestblocklänge kann zunächst bei 20 kb liegen und später anhand von Benchmarkdaten angepasst werden.

### 11.4 Evidenzklassen

```text
high_confidence:
  Coverage über einzigartige Fenster + hochwertige Alignments +
  breakpoint-Evidenz oder mehrere spanning reads

moderate_confidence:
  Stabile Coverage über einzigartige Fenster, aber ohne breakpoint-Evidenz

ambiguous:
  starke Repeat-/Multi-Mapping-Signale, widersprüchliche Evidenz oder
  unzureichende globale Long-Read-Tiefe
```

---

## 12. Phase VII: SV-Calling mit Long Reads

Zusätzlich zur Coverage wird pro Testisolat breakpoint-basiertes SV-Calling durchgeführt.

```bash
sniffles \
  --input results/mapping/TEST01.panel.bam \
  --vcf results/pav/TEST01.sv.vcf.gz \
  --snf results/pav/TEST01.snf \
  --reference results/panel/Mo_multiref_panel_v1.fa \
  --threads 24
```

Für die Kohorte:

```bash
sniffles \
  --input results/pav/TEST01.snf \
          results/pav/TEST02.snf \
          results/pav/TEST03.snf \
          results/pav/TEST04.snf \
          results/pav/TEST05.snf \
          results/pav/TEST06.snf \
          results/pav/TEST07.snf \
          results/pav/TEST08.snf \
          results/pav/TEST09.snf \
          results/pav/TEST10.snf \
  --vcf results/pav/pilot10.cohort.sv.vcf.gz \
  --threads 24
```

Sniffles2 kann Long-Read-basierte Deletionen, Insertionen, Duplikationen, Inversionen und Translokationen detektieren und erlaubt die populationsweite Zusammenführung über `.snf`-Dateien. [web:20][web:22]

Die finale PAV-Entscheidung kombiniert:

```text
PAV-Evidenz = Coverage-Breadth + Mapping-Eindeutigkeit + SV-Breakpoints + spanning reads
```

---

## 13. Phase VIII: Kandidatenregionen hinterlegen und zuordnen

### 13.1 Kandidatenregionen vor dem Mapping definieren

Jede Kandidatenregion erhält eine stabile ID:

```text
ACC_0001 ... ACC_n       accessory region
MCHR_0001 ... MCHR_n     Mini-Chromosom oder mChr-Block
STAR_0001 ... STAR_n     Starship-like region
SUBTEL_0001 ... SUBTEL_n subtelomerische dynamische Region
```

Datei: `results/panel/panel_candidate_regions.tsv`

```tsv
candidate_id	panel_id	region_type	panel_start	panel_end	length_bp	representative_reference	reference_isolates	n_reference_genomes	host_groups	core_accessory_class	repeat_fraction	captain_gene_id	cargo_gene_count	annotation_confidence
MCHR_0001	PANEL000104	mini_chromosome	1	1200000	1200000	AG006	AG006;BR62	2	oryza;eleusine	accessory_chromosome	0.58	NA	42	high
STAR_0004	PANEL000105	starship_like	1	250000	250000	BR62	BR62;X12;X14	3	eleusine;wildgrass	starship_like	0.43	BR62_g08765	17	medium
ACC_0017	PANEL000221	accessory_region	1	85000	85000	OryzaRef3	OryzaRef3;AG006	2	oryza	private_accessory	0.31	NA	6	high
```

### 13.2 Ergebnis pro Testisolat

Datei: `results/pav/candidate_region_calls.tsv`

```tsv
test_isolate	candidate_id	region_type	status	confidence	breadth_5x	mean_depth	unique_window_fraction	sv_support	shared_reference_isolates	closest_panel_reference	assignment_confidence	notes
TEST01	MCHR_0001	mini_chromosome	present	high	0.94	31.2	0.89	12	AG006;BR62	AG006	medium	mChrA-like region
TEST01	STAR_0004	starship_like	absent	high	0.03	0.4	0.96	0	BR62;X12;X14	NA	NA	no interior or boundary support
TEST02	ACC_0017	accessory_region	ambiguous	low	0.43	6.1	0.28	1	OryzaRef3;AG006	OryzaRef3	low	repeat-associated multi-mapping
```

Die Spalten `shared_reference_isolates` und `closest_panel_reference` müssen unterschieden werden:

- `shared_reference_isolates`: Alle Referenzisolate, die den homologen Regioncluster tragen.
- `closest_panel_reference`: Die Referenzvariante mit bester Alignment-Evidenz.
- `assignment_confidence`: Sicherheit dieser spezifischen Zuordnung.

Bei nahezu identischen Regionen ist „mChrA-like und in AG006/BR62 geteilt“ belastbarer als die Behauptung eines eindeutigen donorspezifischen Ursprungs.

### 13.3 Lokalisierung und biologische Annotation

```bash
bedtools intersect \
  -a results/pav/TEST01.pav_blocks.bed \
  -b results/panel/panel_regions.bed \
  -wa -wb \
  > results/pav/TEST01.pav_blocks_localized.tsv

bedtools intersect \
  -a results/pav/TEST01.pav_blocks.bed \
  -b results/panel/panel_genes.gff3 \
  -wa -wb \
  > results/pav/TEST01.pav_gene_overlap.tsv
```

Jeder finale PAV-Call soll enthalten:

```text
PAV_ID
Testisolat
Panel-Koordinate
Länge
PAV-Status
Konfidenzklasse
Regionstyp
Core/accessory-Klasse
Subtelomerisch ja/nein
Repeat-Anteil
überlappende Gene
Orthogruppen
Funktionsannotation
Starship-/Captain-Evidenz
Referenzcluster
geteilte Referenzisolate
nächstähnliche Panelreferenz
```

---

## 14. Phase IX: Rarefaction und Panel-Sättigung

Die Rarefaction bewertet zwei unterschiedliche Fragen:

1. Wie viele neue Kandidatenregionen werden mit zusätzlichen Referenzgenomen entdeckt?
2. Deckt das Referenzpanel die Kandidatenregionen der Long-Read-Testisolate ausreichend ab?

### 14.1 Referenz-Rarefaction

Für jede Panelgröße \(k = 1, ..., 14\):

1. Ziehe viele zufällige Kombinationen von \(k\) Referenzgenomen.
2. Vereinige die in der Kombination vorhandenen Kandidatenregioncluster.
3. Zähle nichtredundante Regionen getrennt nach accessory, Mini-Chromosomen und Starship-like Regionen.
4. Berechne Mittelwert, Median, Standardabweichung und 95%-Intervall.

Formell:

\[
R(k) = \left| \bigcup_{i \in S_k} C_i \right|
\]

Dabei ist \(C_i\) die Menge der Kandidatenregionen der Referenz \(i\), \(S_k\) eine Auswahl von \(k\) Referenzen und \(R(k)\) die Zahl akkumuliert entdeckter nichtredundanter Kandidatenregionen.

Vorschlag für `config/thresholds.yaml`:

```yaml
rarefaction:
  n_permutations: 1000
  saturation_delta_fraction: 0.05
```

Erwarteter Output: `results/rarefaction/rarefaction_reference_panel.tsv`

```tsv
panel_size	region_class	mean_n_regions	median_n_regions	ci_lower	ci_upper	n_permutations
1	accessory	24	24	18	30	1000
1	starship_like	3	3	1	5	1000
...
14	accessory	218	218	218	218	1
14	starship_like	31	31	31	31	1
```

### 14.2 Sättigungskriterium

Berechne den zusätzlichen Gewinn der letzten Panelgröße:

\[
\Delta R_{13 \rightarrow 14} =
\frac{R(14)-R(13)}{R(14)}
\]

Ein Panel wird als näherungsweise gesättigt bewertet, wenn:

- der mittlere Zugewinn beim Hinzufügen der letzten Referenzen unter etwa 2–5 % liegt,
- die Konfidenzintervalle eng werden,
- und die Long-Read-Testisolate nur selten neue, längere und gut gestützte Kandidatensequenzen ohne Panelhomologie liefern.

Berichte getrennte Kurven für:

- accessory Regionen insgesamt
- Mini-Chromosomen
- Starship-like Regionen
- subtelomerische dynamische Regionen
- alle Kandidatenregionen kombiniert

Ein Plateau der Gesamtkurve bedeutet nicht automatisch Sättigung für Mini-Chromosomen oder Starship-like Regionen; gerade diese Klassen können wesentlich langsamer akkumulieren.

### 14.3 Testisolat-Rarefaction und Novelty-Check

Für Testisolate mit ausreichender Sequenziertiefe:

1. Extrahiere unmapped oder schlecht gemappte Reads.
2. Führe gezielte lokale Assemblies durch.
3. Suche die resultierenden Contigs gegen das Multi-Referenzpanel.
4. Definiere gut gestützte, panel-externe Sequenzen als neue Kandidatenregionen.

Eine Novel Candidate Region könnte beispielsweise erfüllen:

- lokaler Contig mindestens 10–20 kb
- nicht überwiegend repetitive Sequenz
- keine Panelhomologie über mindestens 80 % der Länge
- accessory-, mini-chromosomen-, Starship-like- oder subtelomerische Kandidateneigenschaft

Erzeuge dann eine zweite Akkumulationskurve:

```text
x-Achse: Anzahl analysierter Testisolate
y-Achse: kumulative Zahl neuer, im Panel nicht repräsentierter Kandidatenregionen
```

Wenn diese Kurve rasch weiter ansteigt, ist das Panel unvollständig. Wenn sie früh abflacht, spricht das für gute Abdeckung des relevanten Kandidatenraums.

---

## 15. Konfigurierbare Parameter

Datei: `config/thresholds.yaml`

```yaml
panel:
  strict_core_fraction: 1.00
  soft_core_fraction: 0.93
  shell_min_fraction: 0.21
  candidate_min_length_bp: 10000
  dedup_identity: 0.98
  dedup_coverage: 0.90

mapping:
  preset_ont: map-ont
  preset_hifi: map-hifi
  min_mapq_unique: 20
  min_primary_alignment_bp: 3000

pav:
  window_size_bp: 10000
  high_resolution_window_size_bp: 2000
  min_depth_present: 5
  present_breadth_core: 0.80
  present_breadth_accessory: 0.70
  absent_breadth: 0.10
  min_block_bp: 20000
  min_global_depth: 20

starship:
  min_region_length_bp: 20000
  duf3435_evalue: 1.0e-5
  require_cargo_genes: true
  require_boundary_evidence_for_high_confidence: true

rarefaction:
  n_permutations: 1000
  saturation_delta_fraction: 0.05
```

Alle Schwellen sind Startwerte und müssen mit Benchmarkdaten kalibriert werden. Änderungen an Parametern erzeugen eine neue Analyseversion und werden in der Provenance-Datei dokumentiert.

---

## 16. Validierung vor Skalierung auf 42 Isolate

Vor der vollständigen Analyse der übrigen 32 Long-Read-Isolate:

1. Nutze mindestens ein bis mehrere Isolate mit bekannter hochwertiger Assembly als Benchmark.
2. Entferne die jeweilige Assembly temporär aus dem Referenzpanel.
3. Mappe die Long Reads dieses Isolats gegen das verbleibende Panel.
4. Vergleiche die PAV-Calls mit Assembly-vs-Assembly-Analysen.
5. Berechne Sensitivität, Präzision und F1 getrennt für core, accessory, mini-chromosomale, subtelomerische und repeat-reiche Bereiche.
6. Prüfe die Stabilität der Befunde bei moderater Variation der Parameter, beispielsweise MAPQ 20 versus 30 und Breadth 0,70 versus 0,80.

Zielwerte können zunächst sein:

- Sensitivität mindestens 90 % für klar definierte PAVs oberhalb einer Mindestgröße
- Spezifität mindestens 95 % in eindeutigen, nichtrepetitiven Regionen
- konservativer Umgang mit repeats und subtelomerischen Bereichen

---

## 17. Kosten- und Rechenstrategie

| Phase | Aufwand | Vorgehen |
|---|---:|---|
| FASTA-QC, Header, Metadaten | niedrig | Für alle 14 Referenzen sofort durchführen |
| Annotation und Repeatmasking | mittel bis hoch | Einmalig und konsistent für alle Referenzen |
| Whole-genome Alignments/SyRI | mittel | Einmal pro sinnvoller Referenzkombination |
| Panel- und Kandidatenbau | mittel | Versioniert; nur bei Panelupdate wiederholen |
| Mapping der 10 Pilotisolaten | mittel bis hoch | Pilot vor kompletter 42er-Kohorte |
| Coverage-basierte PAV | niedrig bis mittel | Für sämtliche Pilotisolate |
| SV-Calling | mittel | Nach erfolgreichem Mapping und QC |
| Lokale Assemblies unmapped Reads | hoch | Nur gezielt bei Novelty-/HGT-/PAV-Kandidaten |
| Skalierung auf 42 Isolate | skalierbar | Erst nach Validierung der Pilotparameter |

---

## 18. Entscheidungslogik für finale Befunde

Eine biologische Aussage über eine Kandidatenregion sollte folgende Kette erfüllen:

```text
Panelregion definiert
→ Regionklasse und Referenzcluster bekannt
→ ausreichende Long-Read-Qualität im Testisolat
→ Coverage über einzigartige Fenster
→ Bewertung von Multi-Mapping und Repeats
→ optional SV-/Breakpoint- oder spanning-read-Evidenz
→ PAV-Status und Konfidenzklasse
→ Zuordnung zu allen teilenden Referenzisolaten
```

### Interpretation

| Befund | Zulässige Interpretation |
|---|---|
| High-confidence Presence eines Mini-Chromosom-Clusters | Das Testisolat trägt eine mChr-like Region, die mit den im Cluster gelisteten Referenzisolaten homolog ist |
| Presence eines Starship-like Clusters einschließlich Boundary-Support | Das Testisolat trägt wahrscheinlich eine homologe Starship-like Region; strukturelle Bestätigung kann zusätzliche lokale Assembly erfordern |
| Absence einer einzigartigen accessory Region bei ausreichender Tiefe | Robuster Verlust bzw. Nichtvorhandensein im Testisolat |
| Niedrige Coverage in TE-reicher, mehrfach homologer DNA | Ambiguität; keine harte Abwesenheitsaussage |
| Bester Treffer auf einem Referenzisolat, aber mehrere Referenzen im Cluster | Nächstähnliche Referenzvariante, aber keine zwingende donor- oder Ursprungszuweisung |

---

## 19. Abschließendes Prinzip

Das Multi-Referenzpanel ist nicht nur eine FASTA-Datei. Es ist ein versionierter Katalog homologer Sequenzregionen mit:

- stabilen Region- und Cluster-IDs,
- Koordinaten im analytischen Panel,
- ursprünglichen Koordinaten in den 14 Referenzgenomen,
- Core-/accessory-/Mini-Chromosom-/Starship-like-Klassifikation,
- Repeat- und Genannotation,
- Prävalenz im Referenzsatz,
- und einer vollständigen Liste aller Referenzisolate, die die Region teilen.

Dadurch kann jeder Long-Read-PAV später nachvollziehbar als core, accessory, mini-chromosomal, subtelomerisch oder Starship-like lokalisiert werden. Gleichzeitig lässt sich transparent dokumentieren, welche Referenzisolate eine im Testisolat gefundene Kandidatenregion ebenfalls tragen und ob das Panel für die Diversität der untersuchten Blast-Linien bereits ausreichend gesättigt ist.
