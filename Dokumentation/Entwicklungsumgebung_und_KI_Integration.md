# Entwicklungsumgebung und KI-Integration für die Barragan-Pipeline

**Projekt:** Promotion – Barragan / *Magnaporthe oryzae*  
**Datum:** 31. August 2026  
**Zweck:** Aufbau einer reproduzierbaren Arbeitsumgebung für Starfish-, Starship-, Mini-Chromosomen- und Isolatvergleichsanalysen.

---

## 1. Ziel des Setups

Für dieses Promotionsprojekt sollte die Analyse nicht aus einzeln kopierten Terminalbefehlen bestehen. Ziel ist ein versioniertes, reproduzierbares Bioinformatik-Projekt, das von einem einzelnen Testisolat auf viele Isolate erweitert werden kann.

Die technische Grundidee lautet:

```text
VS Code als Arbeitsoberfläche
        +
Git für Versionskontrolle
        +
Conda/Mamba für Software-Umgebungen
        +
Snakemake für automatisierte Pipelines
        +
Python/R für eigene Datenaufbereitung und Auswertung
        +
KI als kontrollierte Unterstützung für Code, Tests und Dokumentation
```

Die KI kann bei Erklärungen, Codeentwürfen, Fehlersuche und Dokumentation helfen. Die eigentliche wissenschaftliche Verantwortung bleibt jedoch beim Projekt: Eingabedaten, Versionen, Parameter, Tests und biologische Interpretation müssen nachvollziehbar gespeichert und geprüft werden.

---

## 2. Empfohlene Gesamtarchitektur

```text
Arbeitsrechner / Browser
        │
        ├── Visual Studio Code
        │   ├── Editor für Python, Snakemake, YAML und Markdown
        │   ├── integriertes Linux-Terminal
        │   ├── Git-Ansicht und Commit-Historie
        │   └── KI-Chat bzw. Pair-Programming
        │
        └── Remote-SSH, falls die Rechenumgebung ein Server ist
                    │
                    ▼
Linux-Rechenumgebung: MilleniumFalke
        │
        ├── Projekt-Repository: ~/promotion/barragan/
        ├── Conda-/Mamba-Umgebungen
        ├── Snakemake-Workflow
        ├── Rohdaten außerhalb von Git
        ├── Ergebnisdaten außerhalb von Git
        └── Logs, Benchmarks und Reports
```

Wenn `MilleniumFalke` ein separater Linux-Server ist, wird VS Code lokal auf dem Arbeitsrechner installiert und über SSH verbunden. Wenn VS Code direkt auf `MilleniumFalke` läuft, kann der Projektordner direkt lokal geöffnet werden.

---

## 3. Visual Studio Code

### Warum VS Code?

VS Code ist eine geeignete Arbeitsoberfläche, weil es Editor, Terminal, Git-Unterstützung, Debugging und Erweiterungen für Python, YAML, Markdown und Remote-Entwicklung kombiniert.

Für dieses Projekt werden vor allem folgende Dateitypen bearbeitet:

| Dateityp | Verwendung |
|---|---|
| `.py` | kleine Skripte für Formatierung, Validierung, Tabellen und Auswertung |
| `Snakefile`, `.smk` | Regeln der Snakemake-Pipeline |
| `.yaml`, `.yml` | Parameter, Software-Umgebungen und Referenzpfade |
| `.tsv` | Probenmetadaten, Mappingtabellen, Ergebniszusammenfassungen |
| `.md` | Dokumentation, Entscheidungen, Methodenentwürfe und Laborbuch |
| `.sh` | einzelne Hilfsbefehle, sofern sie sinnvoller als Python sind |

### Empfohlene Erweiterungen

| Erweiterung | Zweck |
|---|---|
| Remote - SSH | Öffnet einen Ordner auf einem Linux-Rechner/Server über SSH |
| Python | Python-Code, Interpreter-Auswahl, Linting und Debugging |
| Pylance | Autovervollständigung und Typprüfung für Python |
| Snakemake | Syntax-Hervorhebung und Unterstützung für Snakemake-Dateien |
| YAML | Syntaxprüfung und Bearbeitung von YAML-Konfigurationen |
| Markdown All in One | Markdown-Vorschau, Gliederung und Formatierung |
| GitLens | Erweiterte Git-Historie und Dateiversionen |
| Error Lens | Zeigt Fehler und Warnungen direkt im Editor an |
| Jupyter | Optional für explorative Notebooks; nicht für finale Pipeline-Schritte |

### Projektordner öffnen

Direkt auf dem Linux-System:

```bash
code ~/promotion/barragan
```

Falls der Befehl `code` noch nicht verfügbar ist, kann in VS Code über `File → Open Folder…` der Ordner

```text
~/promotion/barragan
```

geöffnet werden.

---

## 4. Remote-SSH einrichten

Wenn `MilleniumFalke` per SSH erreicht wird, ist Remote-SSH die komfortabelste Arbeitsweise. Der Editor arbeitet dann auf dem Server: Dateien, Terminalbefehle, Conda-Umgebungen und Erweiterungen verwenden direkt die Linux-Rechenumgebung.

### SSH-Alias anlegen

Auf dem lokalen Arbeitsrechner kann die Datei `~/.ssh/config` einen Alias enthalten:

```sshconfig
Host milleniumfalke
    HostName DEIN_SERVER_ODER_DEINE_IP
    User flori
    IdentityFile ~/.ssh/id_ed25519
    ServerAliveInterval 60
```

Danach funktioniert die Verbindung im Terminal mit:

```bash
ssh milleniumfalke
```

In VS Code:

1. Command Palette mit `Ctrl+Shift+P` öffnen.
2. `Remote-SSH: Connect to Host...` auswählen.
3. Den Host `milleniumfalke` auswählen.
4. Nach erfolgreicher Verbindung `~/promotion/barragan` als Ordner öffnen.

---

## 5. Conda und Mamba

### Rolle von Conda/Mamba

Bioinformatik-Tools besitzen oft viele Abhängigkeiten. Conda und Mamba stellen sicher, dass definierte Toolversionen gemeinsam installiert werden können.

- **Conda** verwaltet Software-Umgebungen.
- **Mamba** ist weitgehend kompatibel mit Conda, löst Abhängigkeiten aber oft schneller.
- **Ein Environment pro Aufgabenklasse** verhindert Konflikte zwischen Programmen.

### Vorhandene und geplante Environments

| Environment | Aufgabe |
|---|---|
| `starfish_env` | bestehende Umgebung für Starfish, MetaEuk, HMMER und zugehörige Datenbanken |
| `barragan-workflow` | Snakemake, Python und Steuerung der Pipeline |
| `alignment` | später z. B. minimap2, samtools, MUMmer/nucmer, Mash/Skani |
| `barragan-python` | Tabellenanalyse, Validierung, eigene Skripte und Abbildungen |

### Verfügbarkeit prüfen

```bash
git --version
conda --version
mamba --version
snakemake --version
python3 --version
```

Falls Mamba fehlt, kann es in der Basisumgebung installiert werden:

```bash
conda install -n base -c conda-forge mamba
```

### Workflow-Environment erstellen

```bash
mamba create -n barragan-workflow \
  -c conda-forge -c bioconda \
  snakemake \
  python=3.12 \
  pandas \
  biopython \
  pyyaml \
  pytest \
  ruff \
  git \
  -y
```

Aktivieren:

```bash
conda activate barragan-workflow
```

Starfish bleibt zunächst in `starfish_env`. Langfristig definiert die Pipeline eigene Environment-Dateien, sodass Snakemake die richtige Umgebung pro Regel bereitstellen kann.

---

## 6. Git-Repository

### Warum Git?

Git speichert die Geschichte des Codes und der Konfigurationen. Damit lässt sich beantworten:

- Wann wurde ein Parameter geändert?
- Welche Version eines Skripts erzeugte ein bestimmtes Ergebnis?
- Warum wurde eine Regel ergänzt oder entfernt?
- Wie kann ein funktionierender Stand wiederhergestellt werden?

Git speichert **Code, Konfiguration und Dokumentation**, aber im Regelfall keine großen FASTA-, FASTQ-, BAM- oder vollständigen Ergebnisdateien.

### Repository initialisieren

```bash
cd ~/promotion/barragan
git init
```

### Empfohlene Struktur

```text
barragan/
├── README.md
├── .gitignore
├── .github/
│   └── copilot-instructions.md
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
│   │   ├── minichromosomes.smk
│   │   ├── comparison.smk
│   │   └── report.smk
│   ├── scripts/
│   │   ├── normalize_fasta_headers.py
│   │   ├── normalize_gff_seqids.py
│   │   ├── validate_fasta_gff_ids.py
│   │   ├── calculate_contig_metrics.py
│   │   ├── parse_starfish_results.py
│   │   ├── build_pav_matrix.py
│   │   └── score_hgt_candidates.py
│   └── envs/
│       ├── starfish.yaml
│       ├── alignment.yaml
│       ├── python.yaml
│       └── qc.yaml
├── input/
│   ├── manifests/
│   └── metadata/
├── results/
│   ├── qc/
│   ├── normalized/
│   ├── starfish/
│   ├── minichromosomes/
│   ├── comparisons/
│   └── reports/
├── logs/
├── benchmarks/
├── docs/
│   ├── workflow.md
│   ├── decisions.md
│   ├── data_dictionary.md
│   ├── methods_draft.md
│   └── ai_usage.md
└── tests/
    ├── data/
    └── test_normalization.py
```

### `.gitignore`

Lege möglichst früh eine `.gitignore` an:

```gitignore
# Große Rohdaten und Sequenzdaten
assemblies/
raw_data/
reads/
*.fastq
*.fastq.gz
*.fq
*.fq.gz
*.bam
*.bai
*.cram
*.crai

# Große Bioinformatikdateien
*.fna
*.fa
*.fasta
*.fas
*.faa
*.gff
*.gff3
*.hmm
*.dmnd
*.mmi

# Automatisch erzeugte Ergebnis- und temporäre Dateien
results/
logs/
benchmarks/
temp/
tmp/
.snakemake/

# Python- und Editor-Artefakte
__pycache__/
*.pyc
.ipynb_checkpoints/
.vscode/settings.json
.DS_Store

# Zugangsdaten und Geheimnisse
.env
*.key
```

### Erster Commit

```bash
git add README.md .gitignore config workflow docs tests .github
git commit -m "Initialize reproducible Starfish and minichromosome workflow"
```

Danach in kleinen, fachlich verständlichen Einheiten committen:

```bash
git status
git add workflow/scripts/normalize_fasta_headers.py tests/test_normalization.py
git commit -m "Add FASTA header normalization with ID mapping"
```

---

## 7. Snakemake als Pipeline-Engine

### Aufgabe von Snakemake

Snakemake beschreibt, welche Ergebnisdateien aus welchen Eingaben entstehen. Es erkennt Abhängigkeiten zwischen Schritten, kann mehrere Isolate parallel rechnen und führt fehlende oder unvollständige Schritte erneut aus.

Eine vereinfachte Pipeline:

```text
Original-FASTA
      ↓
FASTA-Header normalisieren
      ↓
GFF-SeqIDs synchronisieren
      ↓
FASTA/GFF-Kompatibilität validieren
      ↓
Starfish: YR-Kandidaten suchen
      ↓
Starship-Kontext analysieren
      ↓
Contig-Metriken und mChr-Kandidaten bestimmen
      ↓
Isolatvergleich: PAV, Homologie, Core-vs-Element-Kontrast
```

### Klein beginnen

Der erste Pilotworkflow soll nur vier Arbeitsschritte enthalten:

1. FASTA-Header normalisieren.
2. GFF-SeqIDs mit derselben ID-Mappingtabelle synchronisieren.
3. Prüfen, ob jede GFF-SeqID in der passenden FASTA vorkommt.
4. Starfish-YR-Annotation starten.

Erst wenn dieser Ablauf für ein Isolat stabil funktioniert, folgen drei Isolate und danach der gesamte Datensatz.

### Beispiel: zentraler Snakefile

Datei `workflow/Snakefile`:

```python
configfile: "../config/parameters.yaml"

include: "rules/headers.smk"
include: "rules/starfish.smk"

rule all:
    input:
        expand(
            "../results/starfish/{sample}/{sample}_YR.filt.gff",
            sample=config["pilot_samples"],
        )
```

### Beispiel: Parameterdatei

Datei `config/parameters.yaml`:

```yaml
project_name: barragan_hgt
separator: "_"

threads:
  starfish: 6

pilot_samples:
  - GCA_004346965.1

starfish:
  profile: /home/flori/miniforge3/envs/starfish_env/db/YRsuperfams.p1-512.hmm
  proteins: /home/flori/miniforge3/envs/starfish_env/db/YRsuperfamRefs.faa
  hmm_evalue: 0.001
```

Absolute Pfade sind für den ersten Pilotlauf akzeptabel, müssen aber dokumentiert werden. Später sollten Referenzdateien in einer kontrollierten Ressourcenstruktur oder über zentrale Pfadkonfigurationen verwaltet werden.

### Dry-run vor realen Läufen

Bevor ein rechenintensiver Workflow gestartet wird:

```bash
cd ~/promotion/barragan
conda activate barragan-workflow

snakemake \
  --snakefile workflow/Snakefile \
  --cores 6 \
  --use-conda \
  --dry-run \
  --printshellcmds
```

Ein Dry-run zeigt geplante Befehle und fehlende Eingaben, ohne Dateien zu erzeugen.

Ein realer Lauf:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --cores 6 \
  --use-conda \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds
```

| Option | Bedeutung |
|---|---|
| `--cores 6` | maximal sechs Rechenkerne verwenden |
| `--use-conda` | pro Regel die deklarierte Software-Umgebung nutzen |
| `--rerun-incomplete` | unvollständige Resultate sauber neu erzeugen |
| `--keep-going` | unabhängige Isolate weiterrechnen, wenn eines scheitert |
| `--printshellcmds` | ausgeführte Terminalbefehle sichtbar machen |

---

## 8. KI-Integration

### Geeignete KI-Aufgaben

KI ist besonders hilfreich bei:

- Erklärung von Fehlermeldungen, Tool-Hilfe und Dateiformaten
- Entwurf kleiner Python-Skripte
- Umwandlung eines getesteten Terminalbefehls in eine Snakemake-Regel
- Schreiben und Verbessern von Unit-Tests
- Code-Review und Suche nach Sonderfällen
- Erstellen von YAML-, TSV- und Markdown-Vorlagen
- Formulieren von Methodendokumentation
- Erkennen fehlender Input-/Output-/Log-Deklarationen in Workflow-Regeln

### Aufgaben, die immer selbst geprüft werden müssen

KI darf nicht ungeprüft entscheiden über:

- biologische Interpretation einzelner Treffer
- endgültige HGT- oder Starship-Calls
- Schwellenwerte ohne Kalibrierung
- Löschen, Umbenennen oder Überschreiben großer Datenbestände
- Toolparameter ohne lokale Hilfe, Literatur oder Testlauf
- private Rohdaten, Zugangsdaten, API-Keys oder unveröffentlichte sensible Metadaten

### Sichere Arbeitsreihenfolge

```text
KI-Vorschlag
      ↓
Lokale Tool-Hilfe / Primärquelle prüfen
      ↓
Test mit kleinem Datensatz
      ↓
Output, Logs und Sonderfälle validieren
      ↓
Code testen und reviewen
      ↓
Änderung committen
      ↓
Erst dann auf viele Isolate skalieren
```

### GitHub Copilot in VS Code

Bei Nutzung von GitHub Copilot kann eine projektspezifische Regeldatei angelegt werden:

```bash
mkdir -p .github
code .github/copilot-instructions.md
```

Empfohlener Inhalt:

```markdown
# Barragan bioinformatics workflow instructions

## Scientific and reproducibility rules
- Never alter raw input files in place.
- Write derived files only to results/ or temporary directories.
- Every rule must declare explicit inputs, outputs, logs, threads, and conda environment.
- Do not infer biological conclusions from a single YR hit.
- Treat YR hits as YR_candidate unless additional evidence is documented.
- Preserve isolate IDs and record every ID transformation in a mapping table.
- For FASTA header changes, generate old-to-new ID mappings.
- For GFF changes, validate that every sequence ID exists in the corresponding FASTA.
- Do not hard-code sample names in scripts; read them from config or samples.tsv.
- Do not use destructive commands without first proposing a safe alternative.
- Add or update tests for each custom Python transformation.
- Use clear docstrings, type hints, TSV outputs, and deterministic sorting.

## Code style
- Use Python 3.12 syntax.
- Use pathlib instead of shell-specific path construction where possible.
- Use pandas only for tabular data; use Biopython for FASTA parsing.
- Fail clearly with actionable error messages.
- Keep scripts small: one biological/data transformation per script.
```

Diese Datei ist ein Leitfaden für die KI. Sie soll verhindern, dass sie Rohdaten verändert, Proben fest in Code einbaut oder vorläufige YR-Treffer als bestätigte Starships bezeichnet.

### Beispielprompts

**Python-Skript für FASTA-Header:**

```text
Schreibe ein Python-3.12-Skript für workflow/scripts/normalize_fasta_headers.py.
Eingabe: FASTA und isolate_id.
Ausgabe: normalisierte FASTA mit Headern <isolate_id>_<original_id> sowie eine TSV-Mapping-Datei old_id, new_id.
Verändere keine Sequenzen. Brich mit verständlicher Fehlermeldung ab, wenn ein Header leer oder doppelt ist.
Erstelle außerdem pytest-Tests mit einer minimalen Test-FASTA.
```

**Snakemake-Regel:**

```text
Erstelle eine Snakemake-Regel normalize_fasta_headers.
Der Input kommt aus config/samples.tsv, Output ist results/normalized/{sample}.fna und results/normalized/{sample}.id_map.tsv.
Schreibe ein Log nach logs/headers/{sample}.log.
Nutze workflow/envs/python.yaml, deklarierte threads und keine hart codierten Isolat-IDs.
Erkläre danach jede Zeile kurz.
```

**Review eines Datenformats:**

```text
Prüfe diese GFF- und FASTA-ID-Mapping-Logik kritisch.
Liste mögliche Fehlerfälle auf: FASTA-Header mit Leerzeichen, doppelte Contig-IDs, GFF-SeqIDs ohne FASTA-Gegenstück, Kommentare in der GFF und komprimierte Eingabedateien.
Schlage konkrete Tests vor, aber ändere keine Dateien.
```

---

## 9. Dokumentation und digitales Laborbuch

Neben der Pipeline braucht das Projekt ein nachvollziehbares Methoden- und Entscheidungstagebuch.

### `README.md`

Der README enthält:

- Ziel der Pipeline
- Voraussetzungen und Installation
- kurze Ordnerübersicht
- Beispiel für einen Testlauf
- wichtigste Outputs

### `docs/workflow.md`

Dokumentiert:

- Reihenfolge aller Pipeline-Regeln
- Input und Output jeder Regel
- verwendete Tools und deren Zweck
- erwartete Ressourcen und Laufzeiten

### `docs/decisions.md`

Dokumentiert methodische Entscheidungen, z. B.:

```markdown
## 2026-08-31 — Eindeutige Contig-IDs

**Entscheidung:** FASTA-Header werden zu
`<assembly_accession>_<original_contig_id>` normalisiert.

**Begründung:** Starfish erwartet eine Genome-ID plus Feature-/Contig-ID.
Eindeutige IDs verhindern Kollisionen zwischen Isolaten.

**Konsequenz:** Zugehörige GFF-SeqIDs werden über dieselbe Mapping-Tabelle
synchron angepasst und anschließend gegen die FASTA validiert.
```

Zusätzlich dokumentieren:

- Starfish-Version
- Starfish-Datenbankdateien
- HMM-Schwellenwerte
- MetaEuk-Parameter
- Kriterien für mChr-Klassen
- Kriterien für Starship-Konfidenz
- Kriterien für geteilte Elemente und HGT-Priorisierung

### `docs/ai_usage.md`

Hier wird KI-Nutzung transparent festgehalten:

```markdown
# KI-Nutzungsprotokoll

## Grundsatz
KI dient für Erklärung, Codeentwürfe, Dokumentation und Review.
Alle wissenschaftlichen Entscheidungen, Parameter, Tests und Interpretationen
werden von mir geprüft und dokumentiert.

## 2026-08-31
- Verwendung: Erklärung des Starfish-Workflows und Entwurf einer
  Dokumentationsstruktur.
- Verifikation: Lokale `starfish annotate --help`-Ausgabe geprüft;
  Einzelisolat-Test ausgeführt.
- Ergebnis: 13 HMM-validierte YR-Kandidaten in GCA_004346965.1.
- Interpretation: Als YR_candidate dokumentiert, nicht als bestätigtes Starship.
```

### Datenprovenienz und Prüfsummen

Für jede Assembly sollten Accession, Herkunft, Download-Datum, Datei, GFF-Version und Prüfsumme gespeichert werden.

Beispiel:

```bash
sha256sum \
  assemblies/ncbi_dataset/ncbi_dataset/data/GCA_004346965.1/GCA_004346965.1_ASM434696v1_genomic.fna \
  > input/manifests/GCA_004346965.1.sha256
```

---

## 10. Konkreter Startplan

### Heute: Grundstruktur erstellen

```bash
cd ~/promotion/barragan

mkdir -p \
  config \
  workflow/rules \
  workflow/scripts \
  workflow/envs \
  input/manifests \
  input/metadata \
  results \
  logs \
  benchmarks \
  docs \
  tests/data \
  .github

touch \
  README.md \
  .gitignore \
  config/samples.tsv \
  config/parameters.yaml \
  workflow/Snakefile \
  docs/decisions.md \
  docs/ai_usage.md \
  .github/copilot-instructions.md

git init
```

Danach den Projektordner in VS Code öffnen:

```bash
code ~/promotion/barragan
```

### Diese Woche: ein reproduzierbarer Einzelisolat-Pilot

1. `config/samples.tsv` für `GCA_004346965.1` vervollständigen.
2. Python-Skript zur FASTA-Header-Normalisierung schreiben und testen.
3. Mapping-Tabelle `old_id → new_id` erzeugen.
4. Passende NCBI-GFF suchen und deren SeqIDs mit derselben Mapping-Tabelle anpassen.
5. FASTA/GFF-Validierung implementieren.
6. Starfish mit den normalisierten Eingaben wiederholen.
7. Log, Parameter, Ergebnis und Interpretation als `YR_candidate` dokumentieren.
8. Nach jeder funktionierenden Stufe einen Git-Commit erstellen.

### Danach: Drei-Isolat-Pilot

1. Zwei kontrastierende zusätzliche Isolate in `samples.tsv` aufnehmen.
2. Den Workflow ohne manuelle Pfadänderungen ausführen.
3. YR-Kandidaten in einer Gesamttabelle zusammenführen.
4. Contig-Längen, GC-Gehalt und Genanzahlen als erste mChr-Metriken berechnen.
5. Erste paarweise Alignments und PAV-Logik mit den drei Isolaten testen.
6. Erst nach technischer und biologischer Qualitätskontrolle auf den gesamten Datensatz skalieren.

---

## 11. Was vorerst nicht nötig ist

Für den Einstieg nicht erforderlich:

- Docker oder Kubernetes
- eine komplexe Datenbanklösung
- ein umfassender Machine-Learning-Stack
- KI-Agenten mit Schreibrechten auf alle Datenordner
- vollautomatische HGT-Entscheidungen

Zunächst reichen klare TSV/YAML-Dateien, ein kleines Snakemake-Projekt, geprüfte Python-Skripte und konsequente Dokumentation.

---

## 12. Wichtigster erster Meilenstein

Der erste technische Erfolg ist nicht, alle Isolate gleichzeitig zu analysieren. Der Erfolg lautet:

> Ein Isolat läuft automatisiert von der Original-FASTA über normalisierte, validierte IDs bis zu einem nachvollziehbaren Starfish-YR-Output durch Snakemake.

Erst wenn dieser Ablauf stabil, getestet und dokumentiert ist, sollte die Pipeline auf mehrere Isolate und anschließend auf den vollständigen Datensatz erweitert werden.
