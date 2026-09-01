# Methodische Entscheidungen

## 2026-09-01 — Projekt-Reset auf Basis des POC-Workflow-Dokuments

**Entscheidung:** Das gesamte bisherige Pipeline-Setup (Assembly-basierte
Contig-Klassifikation: `normalize_fasta_headers.py`, Starfish-Integration,
Repeat-/Telomer-/Synteny-Heuristiken, Multi-Referenzpanel für
Core-Konsens, BUSCO/QUAST-QC) wurde entfernt. Neuaufbau auf Basis von
`Dokumentation/POC Workflow – HGT Machbarkeitsstudie M. oryzae.md`, das an
das Dissertationsexposé "Population Genomics of Horizontal Gene Transfer
in Magnaporthe oryzae" (WP1.2–WP1.4) gekoppelt ist.

**Begründung:** Das POC-Dokument beschreibt einen grundlegend anderen
methodischen Ansatz als das, was zuvor iterativ aufgebaut wurde:
Kurzread-Mapping gegen ein Multireferenzpanel mit Coverage-/PAV-basierter
Kandidatendetektion (bwa-mem2/mosdepth), ergänzt um Starfish/Stargraph
(strukturbasiert) und PCA/ARI-Clustering (Host-Lineage-Diskordanz als
HGT-Evidenz) — nicht die zuvor gebaute Assembly-interne
Contig-Heuristik. Beide Ansätze parallel im selben Repo zu pflegen hätte
Verwirrung gestiftet (Nutzerwunsch: "es sollen nur relevante
Informationen zum POC hinterlegt sein").

**Was erhalten blieb:**
- Git-Historie (alle vorherigen Commits weiterhin abrufbar)
- `Dokumentation/` (alle Planungsnotizen, inkl. des neuen POC-Dokuments)
- `assemblies/` (4,1 GB bereits heruntergeladene NCBI-Genome, nicht in
  Git, reiner lokaler Cache)
- NCBI-Isolat-/Sequenzlisten in `data/` und `input/metadata/`

**Was neu aufgebaut wurde (siehe unten für Details):**
- Repo-Struktur exakt nach der im POC-Dokument empfohlenen Gliederung
  (`data/references/`, `data/isolates_poc/`, `envs/`, `workflow/`,
  `results/`, `power_analysis/`)
- Vier Environment-Definitionen (`envs/mapping.yaml`, `starfish.yaml`,
  `python.yaml`, `r.yaml`)
- Snakemake-Grundgerüst (`workflow/Snakefile` + `rules/{mapping,pav,
  starfish,cluster}.smk`), exakt an die im Dokument vorgegebenen Befehle
  angelehnt, lauffähig sobald `config/samples.tsv` und
  `data/references/panel_manifest.tsv` befüllt sind
- Vier neue Skripte: `pav_call.py`, `build_pav_matrix.py`,
  `cluster_and_score.py`, `validate_clustering_synthetic.py`
- `power_analysis/glmm_power_sim.R` (aus dem Dokument übernommen)

## 2026-09-01 — Umfassender NCBI-Sequenzdatensatz (SRA, nicht nur Assemblies)

**Entscheidung:** Zusätzlich zur bereits vorhandenen Assembly-Liste
(`data/ncbi_m_oryzae_assemblies.tsv`/`.jsonl`, ~46 Einträge) wurde der
komplette NCBI-SRA-Katalog für *Pyricularia oryzae* abgerufen
(`data/ncbi_m_oryzae_sra_runs.tsv`, via NCBI E-Utilities
esearch/efetch, `rettype=runinfo`): **3.902 Sequenzierläufe**.
Gefilterte Teilmengen:
- `data/ncbi_m_oryzae_sra_wgs_illumina.tsv` — 1.754 WGS-Kurzread-Läufe
  (Illumina) — Kandidatenpool für den POC-Isolat-Subset (Abschnitt 1)
- `data/ncbi_m_oryzae_sra_wgs_longread.tsv` — 193 WGS-Long-Read-Läufe
  (Nanopore/PacBio) — Kandidatenpool für Referenzpanel-Kandidaten

**Begründung:** Der POC-Ansatz (Kurzread-Mapping gegen Multireferenzpanel)
braucht tatsächliche Rohreads pro Isolat, nicht nur fertige Assemblies.
Die bisherige Assembly-Liste (46 Einträge) deckt nur einen Bruchteil der
tatsächlich bei NCBI verfügbaren *M. oryzae*-Sequenzierdaten ab.

**Verteilung nach Bibliotheksstrategie:** 1.971 WGS, 1.204 RNA-Seq, 435
ChIP-Seq, 86 WCS, 55 Sonstige, 43 WGA, 30 ncRNA-Seq, 25 Bisulfite-Seq, 20
RIP-Seq, 7 Amplicon. Nach Plattform: 3.639 Illumina, 164 Oxford Nanopore,
51 PacBio, 19 DNBSEQ, 17 LS454, 12 BGISEQ.

**Konsequenz:** Vor Auswahl des tatsächlichen POC-Isolat-Subsets (4–6
Referenzen + 10–20 Isolate) müssen aus `ncbi_m_oryzae_sra_wgs_illumina.tsv`
gezielt Isolate nach Host-Lineage-Abdeckung ausgewählt werden (Reis,
Weizen, Fingerhirse, Wildgras je nach Verfügbarkeit) — noch nicht
geschehen, `config/samples.tsv` ist aktuell nur eine leere Vorlage.

## 2026-09-01 — PAV-Calling: Breadth aus mosdepth-Thresholds, nicht nur mittlere Tiefe

**Entscheidung:** `pav_call.py` berechnet die Breadth (Anteil der Region
mit Coverage ≥ `present_depth`) aus mosdepth's `--thresholds`-Ausgabe
(`<prefix>.thresholds.bed.gz`), nicht aus der mittleren Tiefe in
`<prefix>.regions.bed.gz` allein.

**Begründung:** Das POC-Dokument nennt beide Schwellen (Breadth UND
Tiefe) als Kriterium ("80–90 % Breadth bei ≥5× Tiefe"), aber sein
Beispielaufruf übergibt `pav_call.py` nur die `regions.bed.gz`-Datei
(reine mittlere Tiefe pro Region) — daraus lässt sich keine echte Breadth
ableiten (eine Region kann hohe mittlere Tiefe durch einen kleinen,
extrem gut abgedeckten Teilbereich haben, obwohl der Großteil der Region
gar nicht abgedeckt ist). Die `mosdepth --by`-Regel wurde daher um
`--thresholds {present_depth}` ergänzt.

**Konsequenz:** `workflow/rules/pav.smk`s `mosdepth_coverage`-Regel
erzeugt jetzt zwei Ausgabedateien (`.regions.bed.gz` für die Tiefe,
`.thresholds.bed.gz` für die Breadth); `pav_call.py` nimmt beide separat
entgegen (`--mosdepth-regions`, `--mosdepth-thresholds`).

## 2026-09-01 — Region-Klassifikation (core/candidate) über eine gemeinsame BED-Datei

**Entscheidung:** `data/references/candidate_regions.bed` bekommt eine
5. Spalte `region_class` (`core` oder `candidate`), statt core-Marker und
Kandidatenregionen in getrennten Dateien/getrennten PAV-Läufen zu führen.
`pav_call.py` reicht diese Spalte durch, `build_pav_matrix.py` splittet
die kombinierte `candidate_table.tsv` danach in zwei Matrizen für
Abschnitt 3 (PCA/ARI-Clustering).

**Begründung:** Das Dokument trennt konzeptionell "Core-Genom-Marker"
(genomweit, außerhalb der Kandidatenregionen) von "Kandidaten-Regionen"
(aus Abschnitt 2), spezifiziert aber keine konkrete Dateistruktur dafür.
Eine gemeinsame BED-Datei mit Klassenspalte vermeidet einen doppelten
Mapping-/mosdepth-Lauf (einmal pro Regionsklasse) und hält die
Zuordnung an einer Stelle nachvollziehbar.

## 2026-09-01 — Synthetischer Clustering-Validierungstest: ehrlicher Befund zur Methodensensitivität

**Entscheidung:** `validate_clustering_synthetic.py` (POC-Dokument
Abschnitt 3, "jetzt schon möglich") läuft als Monte-Carlo-Simulation über
standardmäßig 20 Wiederholungen, nicht als Einzellauf.

**Begründung/Befund:** Ein einzelner Testlauf mit injizierter
lineage-übergreifender Teilung EINER Kandidatenregion (wie im Dokument
wörtlich beschrieben: "15 % der Isolate teilen eine 'fremde'
Kandidatenregion") zeigt je nach Zufalls-Seed stark unterschiedliche
Ergebnisse: bei 5 Kandidatenregionen und 30 % Empfänger-Isolaten zeigten
nur ~30 % der Einzelläufe (Seeds) einen klaren ARI-Abfall (>0,2); der
gemittelte ARI-Abfall über 20 Wiederholungen liegt bei nur ~0,15 — nicht
zuverlässig genug für ein robustes Diagnosekriterium. Bei mehr
Kandidatenregionen (z. B. 10–30) verschwindet das Signal noch stärker,
da unveränderte, streng lineage-treue Regionen die PCA/k-means-Struktur
dominieren.

**Wichtiger Befund für die tatsächliche Methodik (nicht nur ein
Test-Problem):** Ein aggregierter ARI-Wert über die GESAMTE
Kandidaten-Matrix ist nur schwach sensitiv für ein einzelnes echtes
HGT-Signal, wenn es unter mehreren unauffälligen Kandidatenregionen
verdünnt wird. Für die reale Analyse empfiehlt sich daher zusätzlich zur
aggregierten ARI ein **Pro-Region-Diskordanztest** (z. B.
Lineage-Reinheit pro einzelner Kandidatenregion), bevor man sich allein
auf eine gesamthafte PCA/ARI-Kennzahl verlässt. Dies ist noch nicht
implementiert.

**Konsequenz:** Die Machbarkeitsfrage aus Abschnitt 3 des Dokuments
("Pipeline muss injizierte Diskordanz korrekt erkennen") ist mit der
aktuellen, wörtlich im Dokument beschriebenen Methode **nicht robust
erfüllt** — das muss vor einer Förderzusage/Antragstellung transparent
kommuniziert werden, nicht durch günstige Parameterwahl verdeckt werden.

## 2026-09-01 — R-Umgebung für Power-Analyse noch nicht eingerichtet

**Entscheidung:** `power_analysis/glmm_power_sim.R` wurde 1:1 aus dem
POC-Dokument übernommen (lauffähiges Codegerüst), aber **nicht
ausgeführt** — R ist in dieser Umgebung nicht installiert (das
POC-Dokument merkt das selbst explizit an). `envs/r.yaml` definiert die
nötige Umgebung (r-base, r-glmmTMB, r-simr, r-lme4) zur Einrichtung via
`mamba env create -f envs/r.yaml`.
