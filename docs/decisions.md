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

## 2026-09-01 — Drei Bugfixes beim ersten End-to-End-Lauf auf echten Daten

**Bug 1 (Pfad):** `workflow/rules/regions.smk` hatte
`conda: "../envs/python.yaml"` (nur ein `../`), obwohl die Regel-Dateien
unter `workflow/rules/` liegen und die Envs im Projekt-Root
(`envs/python.yaml`) — nötig ist `"../../envs/python.yaml"`, wie in allen
anderen `.smk`-Dateien. Der Lauf brach dadurch bei `core_marker_regions`
mit `WorkflowError: Error recording metadata` ab, nachdem bereits alle 6
Starfish-Läufe fertig waren. Behoben; die anderen drei `.smk`-Dateien
hatten den korrekten Pfad bereits.

**Bug 2 (Ressourcen, OOM):** bwa-mem2 wurde beim Mapping von
`TH3` durch den Linux-OOM-Killer beendet
(`bwa-mem2.avx2` allein: 4,5 GB RSS). Die WSL2-VM war mit ihrem
Default-Limit (50 % des Host-RAM = 7,5 GB von 16 GB) zu knapp bemessen —
derselbe Grundtyp von Problem wie beim früheren BUSCO-OOM. Da der Host
tatsächlich 16 GB RAM hat, wurde `C:\Users\flori\.wslconfig` neu angelegt
(`memory=11GB`, `swap=4GB`) und WSL neu gestartet (`wsl --shutdown`).
Zusätzlich bekamen `map_to_panel`, `map_longread_to_panel` und `fastp_qc`
in `mapping.smk` jetzt einen expliziten `threads:`-Wert (vorher nur
`params.threads`, d. h. Snakemake reservierte pro Job nur 1 Kern und hätte
bei `--cores 6` mehrere 8-Thread-bwa-mem2-Prozesse gleichzeitig zulassen
können — mit `threads:` werden sie stattdessen korrekt serialisiert).

**Bug 3 (fehlende NaN-Behandlung):** `cluster_and_score()` übergab die
PAV-Matrix ungeprüft an `sklearn.decomposition.PCA`, die NaN-Werte
("uncertain"-Calls) grundsätzlich ablehnt. Das fiel in den bisherigen
Tests/synthetischen Validierungen nie auf, weil dort nur sauberes 0/1
verwendet wurde — echte Coverage-Daten enthalten aber routinemäßig
"uncertain"-Calls (mehrdeutige Breadth/Tiefe). Fix: Regionen, die für
JEDE Isolat-Probe NaN sind, werden vor der PCA verworfen (keine
Information); verbleibende NaNs werden mit dem Spaltenmittelwert
(Präsenz-Häufigkeit über die anderen Isolate) imputiert — eine neutrale
Standardbehandlung für fehlende genotypartige Calls. Test ergänzt
(`test_cluster_and_score_imputes_nan_instead_of_crashing`).

## 2026-09-01 — Erster End-to-End-POC-Lauf auf echten Daten: Ergebnisse und Einschränkungen

**Datengrundlage:** 6 Referenzgenome (70-15/Guy11 = Reis, B71 = Weizen,
CD156 = Fingerhirse, US71 = Foxtail-Hirse, M_grisea = Fingerhirsegras als
Schwestertaxon-Negativkontrolle), 5 Kurzread-Isolate (MZ5-1-6, E34 =
Fingerhirse; TH3 = Reis; U167 = Fingerhirsegras; Py13.1.023 = Weizen) + 1
Langread-Isolat (ZM1-2 = Weizen, Nanopore). 92 Starship/YR-Kandidatenregionen
(Starfish `annotate`, alle 6 Referenzen erfolgreich annotiert, 12–18 YR-Loci
pro Genom) + 120 Core-Marker-Fenster (20 pro Referenz) in
`data/references/candidate_regions.bed`.

**PAV-Calling-Ergebnis (`results/pav_calls/candidate_table.tsv`):**

| Region-Klasse | absent | present | uncertain |
|---|---|---|---|
| candidate (Starship) | 415 | 55 | 82 |
| core | 435 | 74 | 211 |

**Clustering (`results/clustering/ari_summary.tsv`):** Aggregierter ARI
für Core- UND Kandidaten-Matrix identisch bei 0,423 (k=4, n=6 Isolate) —
bei nur 6 Isolaten und 4 Klassen ist die Zahl möglicher Partitionen so
klein, dass ein Zusammentreffen wenig aussagekräftig ist; **diese Metrik
ist bei der aktuellen POC-Stichprobengröße nicht belastbar** (siehe auch
die früher dokumentierte generelle Schwäche des aggregierten ARI).

**Pro-Region-ARI (`results/clustering/candidate_per_region_ari.tsv`):** 7
von 30 auswertbaren Starship-Regionen als "discordant" eingestuft
(z=0,5): `starship_B71_10`, `starship_CD156_3`, `starship_CD156_6`,
`starship_CD156_7`, `starship_CD156_12`, `starship_Guy11_3`,
`starship_M_grisea_10`. Bei nur 6 Isolaten ist auch dieser Test statistisch
sehr instabil (siehe frühere Kalibrierungsergebnisse bei n=40) — als
Hinweisliste für gezielte manuelle Nachprüfung geeignet, nicht als
Beleg.

**Cross-Lineage-Präsenz (händische Zusatzauswertung, nicht Teil der
Pipeline):** 13 Fälle, in denen eine Starship-Kandidatenregion aus einem
Referenzgenom einer ANDEREN Wirtslinie als der des Isolats als "present"
markiert wurde (z. B. `E34`, Fingerhirse-Isolat, zeigt Präsenz für
Starship-Loci aus `B71`/Weizen, `US71`/Foxtail-Hirse und
`M_grisea`/Fingerhirsegras) — das wäre das erwartete Bild eines
lineage-übergreifenden mobilen Elements, ist aber angesichts der
Stichprobengröße als Hypothese, nicht als Nachweis zu behandeln. Bei den
Core-Fenstern gibt es mit 34 Fällen sogar noch mehr Cross-Lineage-Präsenz,
konzentriert vor allem auf `Py13.1.023` (21 von 34) — das deutet darauf
hin, dass die aktuellen Core-Marker-Fenster (zufällig verteilte Fenster
auf dem größten Contig, keine geprüften Single-Copy-Orthologen) echte,
biologisch erwartete Konservierung zwischen Referenzen einfangen, nicht
zuverlässig lineage-diagnostisch sind, und außerdem von
Sequenziertiefe-Unterschieden zwischen Isolaten beeinflusst werden.

**Wichtigste Dateneinschränkung:** `U167` (Fingerhirsegras) hat in KEINER
der 212 Regionen (weder core noch candidate) einen "present"-Call — alle
92 Kandidaten- und 120 Core-Calls sind "absent" oder "uncertain", auch für
Regionen aus seiner eigenen Referenzlinie (`M_grisea`). Die Rohdaten für
`U167` (SRR14705972) sind mit 34 + 37 MB komprimierter FASTQ auffällig
klein gegenüber den anderen Isolaten (250 MB – 1,2 GB) — die
Sequenziertiefe reicht für dieses Isolat vermutlich grundsätzlich nicht
aus, um gegen das 6-Genom-Panel verlässliche Presence-Calls zu erzeugen.
`MZ5-1-6` zeigt ein ähnliches, wenn auch weniger extremes Muster (0
"present" bei den Kandidatenregionen). **Für eine belastbare Aussage
braucht es entweder mehr Sequenziertiefe für diese Isolate oder ihren
Ausschluss aus der PAV-basierten Auswertung.**

**Gesamtfazit:** Die Pipeline läuft jetzt vollständig durch (QC → Mapping
→ PAV → Clustering) und liefert reale, interpretierbare Zwischenergebnisse
— aber bei n=6 Isolaten und ungleicher Sequenziertiefe sind weder der
aggregierte ARI noch die Pro-Region-Diskordanzliste als eigenständiger
Beweis für HGT zu werten. Für die Machbarkeitsstudie ist das erwartbar
(POC-Zweck: Pipeline-Funktionsfähigkeit zeigen, nicht bereits
Signifikanz); für einen Förderantrag sollte explizit benannt werden, dass
die Stichprobengröße und Tiefenunterschiede aktuell die statistische
Aussagekraft begrenzen.

## 2026-09-01 — R-Umgebung eingerichtet, Power-Analyse ausgeführt — Ergebnis: 0 % Power in allen drei Szenarien

**Vorgeschichte:** `power_analysis/glmm_power_sim.R` wurde beim
Projekt-Reset 1:1 aus dem POC-Dokument übernommen, aber nicht ausgeführt
(R war nicht installiert). `envs/r.yaml` wurde jetzt via
`conda env create -f envs/r.yaml` eingerichtet (r-base 4.3.3, glmmTMB
1.1.9, simr 1.0.7, lme4 1.1-37) und das Skript zum ersten Mal ausgeführt.

**Bug 1 (Druckausgabe):** `summary(powerSim(...))` liefert ein
`data.frame` mit den Spalten `successes/trials/mean/lower/upper` — **nicht**
`Power`. Das Original-Snippet griff auf `summary(x)$Power` zu, was immer
`NULL` ergab; `sprintf()` mit einem `NULL`-Argument liefert in R
`character(0)`, wodurch die gesamte Ausgabeschleife lautlos NICHTS
druckte (kein Fehler, keine Zeile) — der erste Lauf sah dadurch fälschlich
nach einem sauberen, aber leeren Erfolg aus. Behoben: `s$mean` (× 100 für
Prozent), `s$successes`/`s$trials` mit ausgegeben.

**Bug 2 (methodisch, tiefer):** `powerSim()` wird im Original-Snippet
aufgerufen, ohne den zu testenden Effekt (`effect_rr`) explizit im Modell
zu verankern — es wird einfach der aus EINEM zufällig simulierten
Datensatz geschätzte Koeffizient als Testziel verwendet. `simr` selbst
kennzeichnet das explizit als **"observed power"-Berechnung**
(`observedPowerWarning`) — ein bekannt verzerrtes, hochvarianzbehaftetes
Verfahren, kein echter Simulations-Power-Test für einen SPEZIFIZIERTEN
Effekt. Der eigentlich korrekte `simr`-Workflow (`fixef(fit) <- ...`, um
den Zieleffekt vor der Simulation explizit zu setzen) wurde getestet und
schlägt an einer echten Kompatibilitätslücke fehl: `simr` 1.0.7s
`fixef<-`-Methode dispatcht auf einen S4-Generic, der glmmTMB-Objekte
nicht kennt (`Error in getClass(cl): "glmmTMB" is not a defined class`).
Ein sauberer Fix würde vermutlich einen Wechsel auf
`lme4::glmer`/`glmer.nb` erfordern (im POC-Dokument selbst als Alternative
genannt) — das wurde in dieser Sitzung nicht umgesetzt (Aufwand/Scope).
Nur der Druckfehler wurde behoben; das Ergebnis unten ist die
"observed power"-Näherung, keine rigorose Power-Analyse.

**Ergebnis (200 Simulationen je Szenario):**

| Szenario | Felder×Jahre×Isolate/Feld-Jahr | Power (Donor×Starship-Interaktion) | Erfolgreiche Sims |
|---|---|---|---|
| Konservativ | 5×2×10 (100 Isolate) | **0,0 % (0,0–1,8 %)** | 0/200 |
| Moderat | 10×2×20 (400 Isolate) | **0,0 % (0,0–1,8 %)** | 0/200 |
| Optimistisch | 15×3×30 (1.350 Isolate) | **0,0 % (0,0–1,8 %)** | 0/200 |

Zusätzlich: Im "konservativ"-Szenario wird der Interaktionsterm
`donor:vectorStarship` selbst — also genau der zu testende Effekt — wegen
Rangdefizienz aus dem Modell entfernt (`dropping columns from
rank-deficient conditional model`). Das Design (nur 5×2=10 Feld-Jahr-
Zellen) ist zu dünn besetzt, um die volle Fixed-Effects-Struktur
(`donor*vector + recipient*vector`, 9 Koeffizienten) überhaupt zu
schätzen. Im moderaten/optimistischen Szenario bleibt der Zielterm zwar
im Modell, die Power ist aber trotzdem exakt 0 % über alle 200
Wiederholungen (0 von 200 Simulationen ergaben ein signifikantes
Ergebnis) — deutlich unter dem, was allein durch Zufallsschwankung zu
erwarten wäre.

**Einordnung:** Wegen Bug 2 (observed power statt spezifizierter Effekt)
ist die exakte Zahl "0,0 %" nicht als präzise Powerschätzung zu
interpretieren, sondern als starkes, reproduzierbares Warnsignal: Die
Donor×Vektor-Interaktion ist mit dem aktuellen Modell/Design bei KEINER
der drei POC-Stichprobengrößen zuverlässig nachweisbar — selbst nicht im
"optimistischen" 1.350-Isolate-Szenario. Das deckt sich mit der im
POC-Dokument selbst formulierten Faustregel ("Interaktionseffekte sind
bei so wenigen Clustern typischerweise deutlich schwerer zu erkennen als
Haupteffekte"), fällt hier aber schärfer aus als das Dokument selbst
erwartet hätte.

**Konsequenz / Empfehlung:**
1. Interaktion vorerst nur explorativ behandeln, nicht als konfirmatorische
   Haupthypothese (genau der im Exposé selbst vorgesehene Fallback:
   "stepwise, beginning with simple models").
2. Für eine belastbare Zahl vor Antragstellung: Modell auf
   `lme4::glmer`/`glmer.nb` umstellen (dort funktioniert `simr`s
   `fixef<-` nativ) und die Power-Analyse mit einem explizit
   spezifizierten Zieleffekt wiederholen.
3. Die Feld×Jahr-Zellenzahl ist der eigentliche Flaschenhals (nicht die
   Isolatzahl pro Zelle) — mehr Felder/Jahre bringen hier vermutlich mehr
   als mehr Isolate pro Feld-Jahr. Das spricht für eine Ausweitung der
   räumlich-zeitlichen Diversität der Stichprobe, nicht nur ihrer Größe.
4. Sobald reale HGT-Kandidatenraten aus Abschnitt 1–3 vorliegen,
   `baseline_rate`/`effect_rr` durch beobachtete Werte ersetzen und
   wiederholen — die aktuellen Platzhalterwerte sind Annahmen, keine
   Beobachtungen.

## 2026-09-01 — Globales Screening: Erweiterungspotenzial aus dem SRA-Katalog

**Befund:** Die 15 größten BioProjects im bereits katalogisierten
SRA-Illumina-Datensatz (`data/ncbi_m_oryzae_sra_wgs_illumina.tsv`, 1.754
Läufe) stellen zusammen 71 % aller Läufe. Eine Titel-Stichprobe dieser
Top-15 zeigt, dass mehrere bereits publizierte **populationsgenomische
Multi-Standort-Studien** sind, nicht verstreute Einzeleinsendungen:
SRP592876 (244 Läufe, Fingerhirse-Isolate Ostafrika), SRP132141 (88,
"Population Genomic Analysis of the Rice Blast Fungus... Expansion of
Three Main Clades"), SRP288432 (48, Populationsgenetik Subsahara-Afrika),
ERP109496 (55, "verschiedene Standorte in Afrika"), SRP076116 (43,
Weizenbrand Brasilien/B71-Linie).

**Konsequenz:** Für Abschnitt 3 (PCA/ARI-Screening) lässt sich die
Wirtslinien-Diversität und Isolatzahl direkt aus bereits publizierten
Daten massiv erweitern, ohne neue Feldarbeit. Rechnerisch/technisch kein
Hindernis (n=50 ≈ 5 h Mapping, ≈ 66 GB Speicher, beides auf dem
aktuellen Rechner machbar). Statistischer Referenzpunkt für "ab welcher
Stichprobe realistische Ergebnisse zu erwarten sind": die eigene
synthetische Validierung (`validate_clustering_synthetic.py`, Default
n=40 = 4 Linien × 10 Isolate) erreichte dort 70–95 % Erkennungsrate für
die Pro-Region-Methode — bei n=6 (aktueller POC-Lauf) strukturell nicht
erreichbar. Faustregel: **mindestens ~8–10 Isolate pro Wirtslinie**.

**Wichtige Einschränkung:** Das globale Screening verbessert nur
Abschnitt 3, nicht automatisch die GLMM-Power aus Abschnitt 4 (siehe
oben) — dafür fehlt den meisten SRA-Einträgen die nötige
Feld/Jahr-Struktur. Ob die oben genannten Multi-Standort-Studien genug
Geo-/Zeit-Metadaten in ihren Supplements für eine Annäherung an das
Feld/Jahr-Design liefern, ist ein offener Prüfschritt.
