# Starfish-Setup und erster erfolgreicher YR-Test

**Projekt:** Promotion – Barragan / *Magnaporthe oryzae*  
**Datum:** 31. August 2026  
**Zweck dieses Dokuments:** Reproduzierbare Dokumentation der bisher durchgeführten Schritte zur Vorbereitung und zum ersten Test von Starfish.

> Dieses Dokument unterscheidet bewusst zwischen erfolgreich ausgeführten Schritten, Warnungen und noch offenen Arbeiten. Es dokumentiert den aktuellen Stand; es behauptet noch nicht, dass Starship-Elemente oder Mini-Chromosomen abschließend identifiziert wurden.

---

## 1. Wissenschaftliches Ziel

Das Promotionsprojekt untersucht genomische Variation, horizontale Genübertragung (HGT), Starship-ähnliche mobile Elemente und akzessorische Mini-Chromosomen in Isolaten von *Magnaporthe oryzae*.

Starfish ist ein Werkzeug zur Suche und Annotation großer mobiler Elemente, insbesondere sogenannter **Starships**. Ein zentrales Merkmal vieler Starships ist ein Gen für eine Tyrosin-Rekombinase (YR), das häufig als **Captain-Gen** bezeichnet wird. Ein einzelner YR-Treffer ist jedoch zunächst nur ein Kandidat und kein vollständiger Nachweis eines Starship-Elements.

Die bisherige Analyse hatte deshalb das begrenzte Ziel:

1. Starfish im vorhandenen Conda-Environment testen.
2. Eine NCBI-Assembly als Eingabe vorbereiten.
3. YR-Kandidaten in dieser Assembly per de-novo-Genvorhersage und HMM-Validierung suchen.

---

## 2. Arbeitsumgebung und Verzeichnisse

Verwendetes Conda-Environment:

```bash
starfish_env
```

Wichtige Arbeitsverzeichnisse:

```text
~/promotion/barragan/starfish
~/promotion/barragan/data
~/promotion/barragan/assemblies/ncbi_dataset/ncbi_dataset/data
```

Verwendete Testassembly:

```text
GCA_004346965.1
Assembly name: ASM434696v1
Genome FASTA:
~/promotion/barragan/assemblies/ncbi_dataset/ncbi_dataset/data/GCA_004346965.1/GCA_004346965.1_ASM434696v1_genomic.fna
```

---

## 3. Erzeugung des NCBI-Assembly-Reports

Zunächst wurde aus dem NCBI-JSONL-Assembly-Report eine tabulatorgetrennte Tabelle erstellt:

```bash
cd ~/promotion/barragan/data

./dataformat tsv genome \
  --inputfile "$HOME/promotion/barragan/assemblies/ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl" \
  --fields accession,assminfo-biosample-accession,assminfo-name \
  > "$HOME/promotion/barragan/starfish/input/assembly_data_report.tsv"
```

Die ersten Zeilen sahen aus wie:

```text
Assembly Accession	Assembly BioSample Accession	Assembly Name
GCA_000292605.2	SAMN02981399	PoP131
GCA_003016745.2	SAMN06050113	ASM301674v2
GCA_004346965.1	SAMN10491321	ASM434696v1
```

### Ergebnis und Problem

Der Versuch,

```bash
starfish format-ncbi \
  --report input/assembly_data_report.tsv \
  --assemblies input/test_assemblies.txt
```

zu verwenden, schlug fehl. Starfish konnte bereits die Kopfzeile nicht als Datenzeile parsen.

Es wurde danach eine headerlose Datei mit zwei Spalten erzeugt:

```text
GCA_004346965.1	ASM434696v1
```

Die Tabulatoren wurden mit folgendem Test bestätigt:

```bash
grep '^GCA_004346965.1' input/assembly_data_report_starfish.tsv | cat -A
```

Ausgabe:

```text
GCA_004346965.1^IASM434696v1$
```

`^I` bedeutet: echter Tabulator.

Trotz korrekter Tabulatoren lehnte `starfish format-ncbi` die Datei ab. Daraus folgt:

- Das Problem war nicht die Tabulatortrennung.
- `format-ncbi` erwartet vermutlich eine spezifischere Reportstruktur als die erzeugte zweispaltige Tabelle.
- Dieser Subcommand wurde für den ersten funktionalen Starfish-Test nicht weiterverfolgt.

Die Datei kann erhalten bleiben, ist für den unten dokumentierten direkten `annotate`-Weg aber nicht erforderlich.

---

## 4. Assembly-Datei für Starfish vorbereiten

Die Hilfe von `starfish annotate` zeigt, dass `--assembly` eine zweispaltige TSV-Datei erwartet:

```text
genomeID<TAB>Pfad-zur-Assembly-FASTA
```

Daher wurde folgende Datei erstellt:

```bash
cd ~/promotion/barragan/starfish

printf '%s\t%s\n' \
"GCA_004346965.1" \
"$HOME/promotion/barragan/assemblies/ncbi_dataset/ncbi_dataset/data/GCA_004346965.1/GCA_004346965.1_ASM434696v1_genomic.fna" \
> input/test_assemblies_2col.tsv
```

Inhalt, geprüft mit sichtbaren Steuerzeichen:

```bash
cat input/test_assemblies_2col.tsv | cat -A
```

Ausgabe:

```text
GCA_004346965.1^I/home/flori/promotion/barragan/assemblies/ncbi_dataset/ncbi_dataset/data/GCA_004346965.1/GCA_004346965.1_ASM434696v1_genomic.fna$
```

Diese Datei ist korrekt formatiert: Die erste Spalte enthält die stabile Genome-ID, die zweite den vollständigen FASTA-Pfad.

---

## 5. Vorhandene Starfish-Datenbanken

Im Conda-Environment wurden folgende Datenbankdateien gefunden:

```bash
find "$CONDA_PREFIX/db" -type f | sort
```

```text
YRsuperfamRefs.faa
YRsuperfams.p1-512.hmm
duf3723.hmm
duf3723.mycoDB.faa
fre.hmm
fre.mycoDB.faa
myb.SRG.fa
myb.hmm
nlr.hmm
nlr.mycoDB.faa
plp.hmm
plp.mycoDB.faa
```

Für den ersten Starship-orientierten Suchlauf wurden verwendet:

```text
YRsuperfamRefs.faa          Referenz-Aminosäuresequenzen für YR-Gene
YRsuperfams.p1-512.hmm      HMM-Profile zur Validierung von YR-Superfamilien
```

Ein **HMM** (Hidden Markov Model) ist hier ein statistisches Profil einer Proteinfamilie. Es erkennt konservierte Sequenzmuster, auch wenn die Proteine nicht exakt identisch sind. Die Kombination aus Referenzprotein-Suche und HMM-Validierung reduziert unspezifische Treffer.

---

## 6. Lokale Hilfe für `starfish annotate`

Die lokale Hilfe bestätigte die zentralen Pflichtargumente:

```text
-a, --assembly    2 column tsv: genomeID, path to assembly FASTA
-p, --profile     profile HMM file
-P, --proteins    FASTA file of query amino acid sequences
-x, --prefix      prefix for naming all output files
-i, --idtag       prefix for predicted gene featureIDs
-o, --outdir      output directory
```

Optional kann mit `-g` eine zweispaltige Tabelle aus Genome-ID und GFF-Dateipfad übergeben werden. Das ist später wichtig, um bereits existierende NCBI-Genmodelle mit neu vorhergesagten YR-Genen zusammenzuführen.

---

## 7. Erfolgreicher erster YR-Testlauf

Vor dem Lauf wurden Ausgabe-, temporäre und Log-Verzeichnisse vorbereitet:

```bash
mkdir -p output/test_YR temp/test_YR logs
```

Ausgeführter Befehl:

```bash
starfish annotate \
  --assembly input/test_assemblies_2col.tsv \
  --profile "$CONDA_PREFIX/db/YRsuperfams.p1-512.hmm" \
  --proteins "$CONDA_PREFIX/db/YRsuperfamRefs.faa" \
  --prefix GCA_004346965.1_YR \
  --idtag YR \
  --outdir output/test_YR \
  --tempdir temp/test_YR \
  --threads 6 \
  > logs/annotate_GCA_004346965.1_YR.log 2>&1
```

### Ergebnis

Der Lauf wurde erfolgreich beendet. Relevante Log-Zeilen:

```text
running metaeuk easy-predict for 1 assemblies..
running hmmsearch on metaeuk annotations..
filtering metaeuk annotations based on hmmsearch results..
found 13 new YR genes and 0 YR genes that overlap with 0 existing genes
done
```

**Interpretation:** In der Assembly `GCA_004346965.1` wurden 13 neu vorhergesagte Gene gefunden, deren Proteine das verwendete YR-HMM-Kriterium erfüllten. Diese Gene sind **YR-Kandidaten**. Sie sind mögliche Captain-Gene, aber noch keine bestätigten Starship-Elemente.

### Erzeugte Ergebnisdateien

```text
output/test_YR/GCA_004346965.1_YR.fas
output/test_YR/GCA_004346965.1_YR.filt.fas
output/test_YR/GCA_004346965.1_YR.filt.gff
output/test_YR/GCA_004346965.1_YR.filt.ids
output/test_YR/GCA_004346965.1_YR.filt.old2new.ids
output/test_YR/GCA_004346965.1_YR.gff
output/test_YR/GCA_004346965.1_YR.hmmout
output/test_YR/GCA_004346965.1_YR.hmmout.ids
output/test_YR/GCA_004346965.1_YR.metaeuk.log
```

### Wichtige Dateien

| Datei | Bedeutung |
|---|---|
| `*.gff` | Genomische Koordinaten der MetaEuk-vorhergesagten Gene |
| `*.filt.gff` | Nach HMM-Validierung gefilterte YR-Genmodelle; für die nächsten Schritte wichtiger als die ungefilterte GFF |
| `*.fas` | Aminosäure- oder Sequenz-Ausgabe der vorhergesagten Gene |
| `*.filt.fas` | Sequenzen der gefilterten, HMM-validierten YR-Treffer |
| `*.hmmout` | detaillierte Ausgabe von `hmmsearch` |
| `*.hmmout.ids` | IDs der HMM-Treffer |
| `*.filt.ids` | IDs der final übernommenen YR-Kandidaten |
| `*.metaeuk.log` | detailliertes Log der de-novo-Genvorhersage |

---

## 8. Header-Warnungen

Während des erfolgreichen Laufs erschienen Warnungen wie:

```text
warning: CP034210.1 ... is being parsed into <2 components using separator '_'.
Make sure ALL sequence headers are formatted like <genomeID><separator><featureID>
```

Die ursprünglichen FASTA-Header enthalten offenbar nur die Contig- bzw. Chromosomen-ID:

```text
>CP034210.1
```

Starfish erwartet für Multi-Genome-Workflows hingegen eine eindeutige zusammengesetzte Kennung, zum Beispiel:

```text
>GCA_004346965.1_CP034210.1
```

### Warum das wichtig ist

Bei mehreren Isolaten können dieselben Contig-Namen auftreten. Wenn in jeder Assembly etwa ein Contig `contig_1` oder `CP034210.1` existiert, wäre die Herkunft später nicht mehr eindeutig. Durch das Präfix mit der Genome-ID wird jede Sequenz global eindeutig.

### Konsequenz

Der Testlauf ist gültig als Nachweis, dass Starfish funktioniert. Für die systematische Multi-Isolat-Pipeline sollen jedoch:

1. FASTA-Header standardisiert werden.
2. Zugehörige GFF-Dateien exakt mit denselben SeqIDs angepasst werden.
3. Die YR-Annotation danach erneut ausgeführt werden.

---

## 9. Was bisher belegt ist

- Starfish lässt sich im Conda-Environment `starfish_env` ausführen.
- Die Direktübergabe einer zweispaltigen Assembly-TSV an `starfish annotate` funktioniert.
- MetaEuk und HMMER wurden erfolgreich auf der Testassembly ausgeführt.
- In `GCA_004346965.1` wurden 13 HMM-validierte YR-Kandidaten gefunden.
- Die ursprünglichen NCBI-FASTA-Header sollten für die Multi-Isolat-Analyse normalisiert werden.
- `starfish format-ncbi` wurde noch nicht erfolgreich eingerichtet und ist für den aktuellen direkten Analyseweg nicht erforderlich.

---

## 10. Was noch nicht belegt ist

Folgende Aussagen dürfen aus dem bisherigen Test **nicht** abgeleitet werden:

- Es wurden 13 vollständige Starship-Elemente entdeckt.
- Die 13 YR-Gene liegen tatsächlich in mobilen Elementen.
- Die Kandidaten besitzen Cargo-Gene oder nachvollziehbare Elementgrenzen.
- Ein Kandidat befindet sich auf einem Mini-Chromosom.
- Zwischen Isolaten wurden horizontale Transfers nachgewiesen.

Für diese Aussagen sind die nachfolgenden Workflow-Schritte erforderlich: Standardisierung aller Eingaben, Konsolidierung mit vollständigen Genannotationen, Starship-Kontextanalyse, Mini-Chromosomen-Klassifikation und isolatübergreifende Vergleiche.

---

## 11. Reproduzierbarkeit

Für jeden zukünftigen Lauf sollen mindestens folgende Informationen gespeichert werden:

- Exakter ausgeführter Befehl
- Tool- und Conda-Environment-Versionen
- Input-Datei und Prüfsumme
- Parameter, insbesondere HMM-E-Wert, Threads und MetaEuk-Optionen
- vollständige Logs
- Datum, Analysezweck und Git-Commit des Workflows

Die Logdatei des erfolgreichen Ersttests liegt unter:

```text
logs/annotate_GCA_004346965.1_YR.log
```

Die aktuellen Resultate liegen unter:

```text
output/test_YR/
```

---

## 12. Nächster dokumentierter Arbeitsschritt

Als nächstes wird eine strukturierte, automatisierte Pipeline erstellt. Sie soll:

1. alle FASTA- und GFF-IDs harmonisieren,
2. Starfish pro Isolat reproduzierbar ausführen,
3. Starship-Kandidaten transparent klassifizieren,
4. Mini-Chromosomen unabhängig identifizieren und vergleichen,
5. presence/absence, Sequenzähnlichkeit und phylogenetische Diskordanz zwischen Isolaten auswerten.

Eine detaillierte Anleitung hierfür befindet sich im separaten Dokument **„Starfish_und_MiniChromosomen_Analyseplan.md“**.
