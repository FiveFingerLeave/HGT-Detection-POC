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

## 2026-09-01 — Pro-Region-ARI implementiert und direkt gegen aggregierte Methode verglichen

**Entscheidung:** `cluster_and_score.py` bekommt eine neue Funktion
`per_region_ari()`: statt einer PCA über die gesamte Kandidaten-Matrix
wird für **jede Kandidaten-Region einzeln** der ARI zwischen ihrem
Präsenz/Abwesenheits-Muster und der bekannten Wirtslinie berechnet
(fehlende/"uncertain"-Werte werden pro Region ausgeschlossen, nicht
geschätzt). `validate_clustering_synthetic.py` führt jetzt beide Methoden
(A = aggregiertes PCA/k-means-ARI wie im Dokument beschrieben, B =
Pro-Region-ARI) auf denselben synthetischen Daten aus und vergleicht sie
direkt. Die neue Ausgabe `results/clustering/candidate_per_region_ari.tsv`
steht auch im echten Pipeline-Lauf zur Verfügung
(`workflow/rules/cluster.smk`).

**Ergebnis des direkten Vergleichs (Monte-Carlo, 20 Wiederholungen):**

| Szenario | Betroffene Isolate | Methode A (aggregiert) | Methode B (Pro-Region) |
|---|---|---|---|
| 5 Kandidatenregionen, 30 % Empfänger | 12/40 | 30 % Erkennungsrate | **95 %** Erkennungsrate |
| 30 Kandidatenregionen, 15 % Empfänger (wörtliche Dokument-Parameter) | 6/40 | **0 %** Erkennungsrate | **70 %** Erkennungsrate |

**Wichtiger methodischer Nebenbefund:** Bei 4 Wirtslinien, aber einem nur
binären (0/1) Merkmal pro Region, kann der Pro-Region-ARI strukturell
**nicht** nahe 1,0 liegen — nach dem Schubfachprinzip müssen mindestens
zwei der vier Linien zufällig denselben Bit-Wert teilen, auch ganz ohne
HGT (bei den obigen Läufen: "saubere" Regionen liegen im Schnitt bei ARI
≈ 0,38–0,39). Aussagekräftig ist daher nicht der Absolutwert, sondern der
**relative** Abstand: liegt die tatsächlich manipulierte Region klar unter
dem Durchschnitt der übrigen Regionen? Das ist der Fall (0,17–0,26 vs.
0,38–0,39) und wird über die Erkennungsrate ("korrekt als
diskordanteste Region identifiziert") gemessen, nicht über einen
absoluten Schwellenwert.

**Konsequenz:** Der Pro-Region-Test ist bei beiden getesteten Szenarien
deutlich sensitiver als der aggregierte Ansatz, erreicht aber selbst bei
den Dokument-eigenen Parametern (30 Regionen, nur 6 betroffene Isolate)
nur 70 % Erkennungsrate — ebenfalls kein Selbstläufer. Für den echten
POC-Lauf empfiehlt sich, `candidate_per_region_ari.tsv` als primäres
Diagnosewerkzeug zu verwenden und `ari_summary.tsv` (aggregiert) nur
ergänzend, nicht umgekehrt.

## 2026-09-01 — Pro-Region-ARI operationalisiert: Klassifikation vs. reine Rangfolge

**Entscheidung:** `classify_region_discordance()` in `cluster_and_score.py`
macht aus dem rohen Pro-Region-ARI-Score eine tatsächliche Entscheidung
pro Region (`discordant` / `concordant` / `uncertain`), über einen
robusten modifizierten Z-Score (Median + MAD, nicht Mittelwert + Std —
damit einzelne echte Ausreißer die "Normal"-Baseline nicht selbst
verzerren). Ausgabe: `results/clustering/candidate_per_region_ari.tsv`
bekommt eine `classification`-Spalte; Schwelle konfigurierbar über
`config/parameters.yaml: clustering.discordance_z_threshold`.

**Kalibrierung (Monte-Carlo, Standard-Szenario: 5 Kandidatenregionen, 30 %
Empfänger, 20–30 Wiederholungen):**

| z-Schwelle | Sensitivität | Falsch-Positiv-Rate |
|---|---|---|
| 0,3 | 56,7 % | 2,5 % |
| 0,5 | 56,7–70,0 % | 2,5 % |
| 0,7 | 40,0 % | 0,0 % |
| 1,0 (ursprünglicher Default) | 36,7–45,0 % | 0,0 % |

Default auf **z = 0,5** gesetzt (bestes gefundenes
Sensitivität/Falsch-Positiv-Verhältnis).

**Wichtiger Befund:** Selbst mit kalibriertem Schwellenwert bleibt die
Klassifikation (auf ~57–70 % Sensitivität) deutlich hinter der reinen
Rangfolge-Methode zurück ("ist diese Region die niedrigst-bewertete von
allen?", 95 % Erkennungsrate im selben Szenario). Grund: Bei nur ~5
Kandidatenregionen pro Batch ist die Median-/MAD-Schätzung selbst
statistisch instabil (kleine Stichprobe) — ein Absolut-Schwellenwert kann
bei so wenigen Regionen nicht so zuverlässig kalibriert werden wie ein
einfacher Rang-Vergleich.

**Konsequenz für den echten POC-Lauf:** Bei den für einen POC realistisch
kleinen Kandidatenzahlen (niedriger zweistelliger Bereich, laut
Dokument-Erwartung) ist **Sortierung nach Pro-Region-ARI (aufsteigend)**
das zuverlässigere Diagnosewerkzeug, nicht der feste
Klassifikations-Schwellenwert. Die `classification`-Spalte bleibt als
schnelle Grobfilterung mit sehr niedriger Falsch-Positiv-Rate nützlich,
sollte aber nicht als alleinige Entscheidungsgrundlage dienen.

## 2026-09-01 — R-Umgebung für Power-Analyse noch nicht eingerichtet

**Entscheidung:** `power_analysis/glmm_power_sim.R` wurde 1:1 aus dem
POC-Dokument übernommen (lauffähiges Codegerüst), aber **nicht
ausgeführt** — R ist in dieser Umgebung nicht installiert (das
POC-Dokument merkt das selbst explizit an). `envs/r.yaml` definiert die
nötige Umgebung (r-base, r-glmmTMB, r-simr, r-lme4) zur Einrichtung via
`mamba env create -f envs/r.yaml`.
