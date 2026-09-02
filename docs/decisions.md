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

## 2026-09-01 — Pivot: neues Workflow-Dokument (`POC_HGT_Starships_Workflow.md`), alte Kurzread-Coverage-Pipeline entfernt

**Entscheidung:** Der Nutzer hat ein neues, deutlich umfassenderes
Workflow-Dokument bereitgestellt
(`Dokumentation/POC_HGT_Starships_Workflow.md`, ursprünglich
`C:\Users\flori\Downloads\POC_HGT_Starships_Workflow.md`) und angewiesen,
alles bisher Erstellte zu entfernen, was unter dem neuen Ansatz nicht mehr
gebraucht wird. Das alte Kurzread-Coverage-Pipeline-Setup (`workflow/`,
`config/`, `power_analysis/`, alte `tests/`, `data/references/panel_manifest.tsv`
und `candidate_regions.bed`, die bwa-mem2-Indexdateien) wurde entfernt.

**Was sich inhaltlich ändert:** Das neue Dokument ersetzt den bisherigen
Ansatz (Kurzread-Mapping gegen ein kleines 6-Genom-Panel + Coverage-PAV +
PCA/ARI-Clustering) durch einen 8-Phasen-Workflow, der auf einem größeren
(20–30 Isolate), überwiegend Long-Read-assemblierten und einheitlich
annotierten Referenzpanel aufbaut: Assembly/QC (BUSCO, BlobTools) →
Annotation/Repeat-Masking (funannotate/BRAKER, RepeatModeler/Masker) →
Starship-Katalog (weiterhin `starfish`) → Pangenom-PAV (Panaroo/PIRATE/
Roary + Alignment-basierte akzessorische Regionen) → **Sättigungs-/
Rarefaktionsanalyse als Kernstück des POC** (Subsampling-Kurve: ab welcher
Panelgröße flacht die Kandidatenentdeckung ab?) → Validierung (Leave-one-out
+ echte Kurzread-Testisolate) → phylogenetische Inkongruenz (Core-Baum vs.
Starship-Gen-Baum) → GLMM-Stichproben-/Modellkomplexitäts-Simulation (neu
gerahmt: klonlinienbasierte simulierte Phylogenien statt Feld/Jahr-Design)
→ Synthese/Go-No-Go-Bericht.

**Was erhalten blieb (siehe README.md für Details):**
- Git-Historie (Checkpoint-Commit vor der Bereinigung: alle entfernten
  Dateien weiterhin über `git log`/`git show` abrufbar)
- `Dokumentation/` (alle bisherigen Planungsnotizen + das neue Dokument)
- `assemblies/` (6 bereits heruntergeladene Referenzgenome — Startpunkt für
  das größere 20–30-Genom-Panel aus Phase 1)
- `data/isolates_poc/` (5 Kurzread- + 1 Langread-Testisolat, bereits
  heruntergeladen und QC-geprüft — direkt nutzbar für Phase 5)
- Alle NCBI-Kataloge in `data/` und `input/metadata/`
- `envs/starfish.yaml` (Starfish wird in Phase 2 unverändert weiter
  gebraucht); `envs/mapping.yaml`/`python.yaml` als generische
  Werkzeug-Bausteine (minimap2 wird in Phase 3 für Alignment-basierte
  PAV-Detektion erneut gebraucht)

**Neu angelegt:** `poc_hgt_starships/{00_data...09_report}/` mit
Phasen-READMEs (Input/Output/Tools je Phase, aus dem neuen Dokument
übernommen).

**Bewusst zurückgestellt:** Phase 7 (Stichproben-/GLMM-Simulation) — auf
expliziten Nutzerwunsch ("lass das GLMM erstmal außen vor") nicht jetzt
implementiert, obwohl das neue Dokument sie vorsieht. Die frühere
Analyse zum alten Feld/Jahr-GLMM (siehe Eintrag oben, 0 %-Power-Befund)
bleibt als Hintergrundwissen dokumentiert, ist aber für das neu gerahmte
Phase-7-Modell (klonlinienbasierte Simulation) nicht direkt übertragbar.

**Nächster Schritt:** Phase 1 (`poc_hgt_starships/00_data/`) — Kuration
des 20–30-Isolat-Panels aus den bereits vorhandenen NCBI-Long-Read-
Assembly-Katalogen (`data/ncbi_m_oryzae_longread_candidates*.tsv`,
`ncbi_m_oryzae_longread_highquality.tsv`, `ncbi_m_oryzae_sra_wgs_longread.tsv`),
mit Fokus auf maximale Klonlinien-/Wirtsdiversität.

## 2026-09-02 — Vollständiger NCBI-Assembly-Katalog (605 Genome) via Datasets REST API

**Entscheidung:** Statt der bisherigen ~46-Einträge-Assembly-Liste
(`data/ncbi_m_oryzae_assemblies.tsv`, älterer, unvollständiger Pull) wurde
ein vollständiger Neuabruf über die NCBI Datasets REST API v2
(`https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/318829/dataset_report`,
txid 318829 deckt sowohl "Pyricularia oryzae" als auch das taxonomische
Synonym "Magnaporthe oryzae" ab) durchgeführt: **605 Assemblies**, in
einem einzigen Request (page_size=1000) vollständig abgerufen. Ergebnis in
`poc_hgt_starships/00_data/ncbi_pyricularia_oryzae_assemblies_full.tsv`
(27 Spalten: Accession, Organismus, Stamm, Assembly-Level, `sequencing_tech`
— direkt aus der API, keine Heuristik nötig —, BioSample-Host/-Isolationsquelle/
-Geo/-Sammeldatum, Contig-/Scaffold-N50/-L50, Chromosomenzahl, GC%).

**Warum die Datasets-API statt E-Utilities:** Die Datasets-API liefert
`assembly_stats` (N50 etc.) UND `assembly_info.sequencing_tech` UND
eingebettete BioSample-Attribute in einem einzigen strukturierten
JSON-Response pro Assembly — bei E-Utilities (esummary) hätte das
mindestens 3 getrennte Abfragen pro Assembly gebraucht (Assembly-Summary,
BioSample-Summary, ggf. Assembly-Stats-Report-Datei).

**Kernbefunde (Details: `poc_hgt_starships/00_data/dataset_summary.md`):**
- **Speicher:** 24,2 GB unkomprimierte FASTA-Summe (605 Genome à
  35–49 Mb, Median 39,9 Mb); ≈ 6,1 GB gzip-komprimiert.
- **Host-Diversität:** 29 normalisierte Kategorien nach
  Tippfehler-/Synonym-Bereinigung (`host_diversity_summary.tsv`);
  dominiert von *Oryza sativa* (275), *Triticum aestivum* (59),
  *Eleusine* spp. (41), *Urochloa* spp. (19), *Lolium* spp. (16);
  **27 % (164/605) ohne Host-Attribut im BioSample-Datensatz.**
- **Assembly-Level/Contig-Größe:** nur 14 "Complete Genome" (2,3 %) + 42
  "Chromosome"-Level (6,9 %) — 66 % (398) nur Scaffold-Level mit
  Contig-N50-Median von nur 0,06 Mb.
- **Kein explizit T2T-geflaggtes Assembly** (0 von 605). 14 Assemblies
  sind "gapless chromosome-level" (Contig-Zahl = Chromosomenzahl) — ein
  Näherungskriterium für T2T-Qualität, aber ohne verifizierte
  Telomer-Repeats an beiden Enden, daher nicht mit echtem T2T
  gleichzusetzen.
- **Datentyp:** 84 % Short-Read-basiert (508), 11,6 % Long-Read (70),
  4 % Hybrid (24), Rest Sanger/historisch. Erwarteter Zusammenhang
  bestätigt: Complete-Genome-Level ausschließlich Long-Read/Hybrid,
  Scaffold-Level fast ausschließlich Short-Read.

**Konsequenz für Phase 1:** Von 605 Assemblies kommen nur die 56
Complete-Genome-/Chromosome-Level-Einträge realistisch als
Referenzpanel-Kandidaten in Frage (Rest zu fragmentiert für
Starship-Boundary-Calling, siehe Dokument-Vorgabe "Repeat-Masking
essenziell für saubere Boundary-Calls" — auf Scaffold-Level mit
Contig-N50 ~60 kb sind Starship-Grenzen kaum sauber bestimmbar). Von
diesen 56 sind wiederum nur die mit Host-Attribut UND ausreichender
Klonlinien-Diversität für die 20–30-Isolat-Zielgröße relevant — noch zu
filtern.

## 2026-09-02 — Zweiter Pivot: `multireferenzpanel_pav_workflow.md`, eigenes Projekt `magnaporthe_multiref_pav/`

**Entscheidung:** Der Nutzer hat ein drittes, noch spezifischeres
Workflow-Dokument bereitgestellt
(`Dokumentation/multireferenzpanel_pav_workflow.md`, 19 Abschnitte) und
angewiesen, dessen Arbeitsschritte strikt zu befolgen und dafür eine
eigene Ordnerstruktur zu verwenden. Anders als beim ersten Pivot (siehe
oben) wurde NICHTS aus `poc_hgt_starships/` entfernt — das neue Dokument
ist eine sehr viel konkretere Ausarbeitung speziell des
Multireferenzpanel-/PAV-Teils (deckungsgleich mit
`poc_hgt_starships/{01_assembly_qc,03_starship_calls,04_pav_matrix}`),
nicht ein Ersatz für den gesamten 8-Phasen-Plan. Neues, eigenständiges
Projektverzeichnis `magnaporthe_multiref_pav/` exakt nach der in
Abschnitt 4 des Dokuments vorgegebenen Struktur angelegt
(`config/`, `data/{references,annotations,longreads,resources}/`,
`envs/`, `workflow/{rules,scripts}/`, `results/`, `logs/`).

**14-Genom-Referenzpanel identifiziert und heruntergeladen:** Die im
Dokument geforderten "14 vollständigen Genomassemblies" entsprechen
exakt den 14 "Complete Genome"-Einträgen aus dem tags zuvor erhobenen
605-Genome-NCBI-Katalog (nicht neu gesucht, direkt weiterverwendet).
Alle 14 FASTA erfolgreich über die NCBI Datasets API heruntergeladen
(`data/references_raw/`) — **nur 1 von 14 (`GCA004346965_1`,
Eleusine-Isolat) hat eine mitgelieferte GFF3-Annotation**, die übrigen
13 brauchen die in Abschnitt 6.2 vorgesehene Reannotation (BRAKER3/
Liftoff), bevor genbasierte Analysen (OrthoFinder, Abschnitt 7.1)
möglich sind.

**Kritischer Befund zur "42 Long-Read-Testisolate"-Annahme:** Das
Dokument geht davon aus, dass 42 chromosomenbasierte Long-Read-Isolate
für die Pilotauswahl (Abschnitt 9) zur Verfügung stehen. Der NCBI-Katalog
liefert exakt 42 Assemblies auf Chromosome-Level (Zahlen-Übereinstimmung
kein Zufall) — aber bei genauerer Prüfung der `sequencing_tech`-Metadaten:
- **Nur 15 von 42 sind tatsächlich Long-Read/Hybrid-sequenziert**
  (13 long-read + 2 hybrid); 25 sind Short-Read-basiert (überwiegend eine
  große Charge von 22 brasilianischen Weizen-Isolaten, vermutlich
  referenzgestützt gescaffoldet, nicht de-novo Long-Read-assembliert);
  2 sind das historische Sanger-70-15-Duplikat (GCA/GCF_000002495.2).
- **Host-Diversität in diesem 42er-Pool ist stark verzerrt:** 31/42
  *Triticum*, 5 *Oryza*, je 1 *Lolium*/*Setaria*, **0 *Eleusine*.** Die
  in Abschnitt 9.1 geforderte Stratifizierung (u. a. "2
  Eleusine-assoziierte Isolate") ist aus diesem Pool nicht erfüllbar.

**Konsequenz:** Für echte Eleusine-Long-Read-Testisolate muss der
separate SRA-Rohdaten-Katalog (`data/ncbi_m_oryzae_sra_wgs_longread.tsv`,
193 Läufe) nach BioSample-Host-Attributen durchsucht werden — noch nicht
geschehen. `config/samples_candidate_pool.tsv` (42 Zeilen, mit
`platform`-Spalte) dokumentiert den vollen Pool inkl. dieser
Einschränkung; die finalen 10 Pilotisolate (`config/samples.tsv`) sind
deshalb noch nicht befüllt.

**Phase I (Abschnitt 6.1) ausgeführt:** `seqkit stats` für alle 14
Referenzen (`results/qc/assembly_stats.tsv`) — konsistentes Bild (7–10
Contigs, 42–48 Mb, N50 5,7–7,5 Mb, GC ~50 %, passend zum erwarteten
*M. oryzae*-Profil). BUSCO (`sordariomycetes_odb10`) lief im Hintergrund
für alle 14 Genome (env `qc_env`, bereits aus einer früheren Sitzung mit
BUSCO vorhanden, nur `seqkit` ergänzt statt einer redundanten neuen
Umgebung).

**Bewusst zurückgestellt (Umfang/Werkzeugverfügbarkeit):** Annotation
(BRAKER3 — braucht separate GeneMark-Lizenz, kein reiner Conda-Install),
Repeat-Masking (RepeatModeler2/EDTA), Whole-genome-Alignment/SyRI,
OrthoFinder, Panel-Bau, Long-Read-Mapping, fensterbasierte PAV,
SV-Calling (Sniffles2), Rarefaction. Alle als dokumentierte Stubs in
`workflow/rules/*.smk` mit Status-Kommentar und Voraussetzungen angelegt,
damit die Snakemake-Struktur vollständig ist und der nächste
Implementierungsschritt pro Datei klar ist.

**Bug: BUSCO-OOM bei paralleler Ausführung.** Der erste BUSCO-Lauf
(`--cores 2`, 2 Genome parallel) führte zu wiederholten, zunächst
verwirrenden Abbrüchen (leere `/tmp`, `LockException`, EXT4-Un-/Remount
und "journal corrupted or uncleanly shut down" im Kernel-Log — sah nach
einem WSL-VM-Neustart aus). `dmesg` zeigte die eigentliche Ursache klar:
`Out of memory: Killed process ... (python3) ... anon-rss:6897092kB` —
ein einzelner BUSCO-Genome-Mode-Lauf (metaeuk-Genvorhersage gegen ein
~44-Mb-Genom) braucht allein **~6,9 GB RSS**; zwei parallel sprengen die
10 GB WSL-RAM-Grenze (+4 GB Swap) klar. Fix: `qc.smk`s
`busco_reference`-Regel bekommt jetzt eine explizite `threads:`-Direktive
(reserviert das volle Kernbudget pro Job), und der Lauf wird mit
`--cores 1` (echte Serialisierung, ein Genom nach dem anderen) statt
`--cores 2` gestartet. Erwartete Laufzeit dadurch länger (~14 × 10–20 Min
statt parallelisiert), aber stabil.

## 2026-09-02 — BUSCO-Abbrüche: tatsächliche Root Cause war WSL2 `autoMemoryReclaim`, nicht Energiesparmodus

**Vorgeschichte:** Trotz der OOM-Fixes (siehe oben) brach der serialisierte
BUSCO-Lauf (`--cores 1`) weiterhin unvermittelt ab — ohne internen
BUSCO-Fehler (Abbruch mitten in `hmmsearch`, teils sogar mit einem
WSL-Interop-Fehler "Failed to start the systemd user session"). `dmesg`
zeigte durchgehend EXT4-Un-/Remount-Zyklen der Root-Disk (`sdd`) im
~100–130-Sekunden-Takt, begleitet von
`systemd-journald: File ... corrupted or uncleanly shut down`.

**Erste (falsche) Hypothese:** Windows-Energiesparmodus (AC-Standby-Timeout
45 Min) bzw. USB Selective Suspend (power-cycelt vermeintlich die Disk
hinter `sdd`). Nutzer wurde befragt und wählte "Sleep-Timeout temporär
deaktivieren"; beides (`STANDBYIDLE` und USB Selective Suspend, AC-Seite)
wurde per `powercfg` deaktiviert. **Ergebnis: kein Effekt** — der
identische Abbruch-Rhythmus trat unverändert erneut auf, diesmal mit dem
zusätzlichen Fund, dass `/dev/sdd` gar kein externes/USB-Laufwerk ist,
sondern die **Root-Disk der WSL2-VM selbst** (`/` und
`/mnt/wslg/distro`, 1 TB dynamisches VHDX) — die powercfg-Hypothese war
damit strukturell unplausibel (kein USB-Gerät betroffen) und wurde
verworfen.

**Tatsächliche Root Cause:** `C:\Users\flori\.wslconfig` hatte keine
explizite `autoMemoryReclaim`-Einstellung, wodurch WSL2 (Version 2.7.12.0)
den Default **`gradual`** verwendet — eine periodische
Arbeitsspeicher-Kompaktierung der VM, die bei speicherintensiven
Workloads (wie dem ~6,9 GB RSS BUSCO-Prozess) die VM kurz genug
einfriert, um I/O-Timeouts auf der Root-Disk und dadurch die beobachteten
Remounts/Journal-Korruption auszulösen. Fix: `.wslconfig` um
```
[experimental]
autoMemoryReclaim=disabled
```
ergänzt, `wsl --shutdown` ausgeführt (sauberer Neustart, `uptime` bestätigt
0 Min), BUSCO-Lauf erneut gestartet. **Ergebnis: alle 14 Genome liefen ohne
Unterbrechung durch.** Die powercfg-Änderungen wurden auf die
ursprünglichen Werte zurückgesetzt (AC-Standby 0x00000a8c/2700s,
USB Selective Suspend AC 0x00000001/aktiviert), da sie nachweislich nicht
die Ursache waren.

**Lektion:** Bei WSL2-VM-internen Stabilitätsproblemen (Un-/Remounts der
Root-Disk, nicht eines Peripheriegeräts) zuerst `.wslconfig`
(`autoMemoryReclaim`, `vmIdleTimeout`, `sparseVhd`) prüfen, bevor
Windows-Host-Energieeinstellungen als Ursache vermutet werden — das
Gerät hinter der scheinbar "unmount/remount"-betroffenen Disk sollte
immer zuerst per `mount`/`lsblk` identifiziert werden (hier: `sdd` = `/`,
kein USB-Gerät).

**BUSCO-Ergebnis (`sordariomycetes_odb10`, alle 14 Referenzgenome):**
durchweg 97,9–98,2 % Complete (größtenteils Single-Copy, Duplication
≤0,5 %), <2 % Missing — konsistent hohe Assembly-Vollständigkeit über das
gesamte Panel, keine Ausreißer.

## 2026-09-02 — Eleusine-Long-Read-Lücke geschlossen: 2 echte Isolate im SRA-Rohdatenkatalog gefunden

**Vorgehen:** Aus `data/ncbi_m_oryzae_sra_wgs_longread.tsv` (193 Läufe) die
158 eindeutigen BioSample-Accessions extrahiert und per NCBI E-Utilities
(`efetch db=biosample`, Batches à 50) die BioSample-Attribute (Host,
Isolate, Isolation Source, Geo) abgerufen (`config/longread_sra_
biosample_hosts.tsv` im neuen Projekt gesichert, 115/158 BioSamples
lieferten Attribute — der Rest sind ältere BioSamples ohne strukturierte
Metadaten).

**Ergebnis:** Zwei Isolate mit Host `Eleucine coracana` (Fingerhirse) UND
echten long-read-Rohdaten gefunden:

| Isolat | BioSample | SRA-Run | Plattform | Bases | Herkunft |
|---|---|---|---|---|---|
| K23/123 | SAMN08033374 | SRR6307184 | PacBio RS II | 3,31 Gb (≈74× auf 44,5 Mb) | Kenia: Busia district, Halsbrand |
| E34 | SAMN12142210 | SRR9972918 | PacBio Sequel | 8,41 Gb (≈189×) | Äthiopien: Diga, Halsbrand |

Beide als neue Zeilen an `config/samples_candidate_pool.tsv` angehängt
(Pool jetzt 44 statt 42 Einträge; `fastq`-Spalte trägt bereits den
SRA-Run-Accession, da die Rohdaten anders als bei allen 42 ursprünglichen
Chromosome-Level-Einträgen tatsächlich direkt herunterladbar sind).
`read_n50_bp` ist hier die mittlere Subread-Länge aus der SRA-Runinfo
(keine echte N50, da diese Kennzahl im runinfo-Format fehlt) — vor dem
tatsächlichen Mapping-Lauf sollte die echte N50 aus den heruntergeladenen
Reads berechnet werden.

**Wichtiger Nebenbefund:** Cross-Referenz der 42 ursprünglichen
Chromosome-Level-BioSamples gegen die 158 long-read-SRA-BioSamples ergab
**null Überschneidungen** — keines der 42 Assemblies hat in diesem
Katalog auffindbare Rohreads (die Long-Read-Assemblies wurden offenbar
ohne Rohdaten-Deposit eingereicht, oder unter einem anderen BioSample als
dem Assembly-BioSample). Das heißt: `config/samples_candidate_pool.tsv`
enthielt bislang bei KEINEM der 42 Einträge tatsächlich ladbare FASTQ
("assembly only" bei allen) — K23/123 und E34 sind damit nicht nur der
Eleusine-Fix, sondern aktuell die EINZIGEN beiden Einträge im gesamten
Pool mit real verfügbaren Rohdaten. Für die übrigen Host-Gruppen
(Triticum, Oryza, Lolium, Setaria) muss vor der finalen 10-Pilotisolate-
Auswahl (Abschnitt 9.1) ebenfalls im SRA-Katalog nach Rohdaten gesucht
werden, nicht nur nach Assemblies — noch nicht geschehen.

## 2026-09-02 — Repeat-Masking (Abschnitt 6.3): eigene, schlanke Environment statt der schweren `annotation.yaml`

**Entscheidung:** `envs/repeats.yaml` (neu, `multiref-repeats`) mit nur
`repeatmodeler=2.0.5`, `repeatmasker=4.1.7`, `bedtools`, `seqkit`,
`samtools` angelegt, statt RepeatModeler/RepeatMasker aus der bereits
geplanten `envs/annotation.yaml` zu installieren (die zusätzlich BRAKER3,
liftoff, eggnog-mapper, orthofinder, diamond, iqtree bündelt — ein
gemeinsamer Install-Versuch aller dieser Pakete wäre langsamer und
fragiler gewesen, insbesondere wegen BRAKER3s GeneMark-Lizenzabhängigkeit,
die ohnehin separat behandelt werden muss).

**`workflow/rules/repeats.smk` implementiert** (vorher reiner Stub):
`index_reference_fai` (samtools faidx) → `build_repeat_database`
(BuildDatabase, NCBI-Engine) → `run_repeatmodeler` (mit `-LTRStruct` für
LTR-Retrotransposon-Sensitivität) → `run_repeatmasker` (mit der
genomspezifischen RepeatModeler-Bibliothek) → `repeat_windows`
(`workflow/scripts/repeat_windows.sh`: RepeatMasker-`.out` → BED →
`bedtools merge`/`coverage` gegen ein Fenstergitter). Fenstergröße wird
aus `config/thresholds.yaml: pav.window_size_bp` (10 kb) übernommen,
nicht neu definiert — damit lässt sich die Repeat-Dichte pro Fenster
später direkt mit der windowbasierten PAV-Klassifikation (Abschnitt 8)
joinen, exakt wie im Dokument gefordert ("Repeat-Anteil" als Spalte der
finalen Panel-Tabelle, Abschnitt 8.4/8.5). Dafür musste
`workflow/Snakefile` erweitert werden, um `config/thresholds.yaml`
selbst einzulesen (`thresholds = yaml.safe_load(...)`) — vorher wurde die
Datei nur referenziert, nie geparst.

**Wie bei BUSCO:** `run_repeatmodeler`/`run_repeatmasker` bekommen ein
explizites `threads:` (= `config["threads_default"]`), um bei `--cores`
gleich dem Threadwert eine echte Serialisierung zu erzwingen — RepeatModeler
ist ähnlich speicher-/zeitintensiv wie BUSCO, mehrere parallele Läufe auf
der 10-GB-WSL2-VM sind ein bekanntes Risiko (siehe BUSCO-OOM oben).

**Pilotlauf vor Vollausführung:** Da RepeatModeler2 mit `-LTRStruct` für
Genome dieser Größe laut Literatur mehrere Stunden pro Genom brauchen
kann und 14 Genome seriell potenziell 1–4+ Tage bedeuten, wurde zunächst
NUR `GCA036493215_1` (kleinstes Genom im Panel, 42,5 Mb) als Zeitpilot
gestartet, bevor alle 14 als langer Hintergrundlauf committed werden.

**Pilotlauf-Ergebnis (Zeitmessung):** Runde 1 (RepeatScout, größte
Stichprobe) brauchte 63 Minuten für `GCA036493215_1` — davon allein 56
Minuten für EINE einzelne, ungewöhnlich kopienreiche Repeat-Familie
("family-0"; die übrigen 88 entdeckten Familien liefen in Sekunden bis
niedrigen Minuten durch). RepeatModeler ist dabei, anders als BUSCO,
**kaum speicherhungrig** (~1 GB RSS statt ~7 GB) — die 10-GB-WSL2-Grenze
ist hier kein Thema, der Engpass ist rein CPU-/Zeit-gebunden. Damit ist
echte Parallelisierung über mehrere Genome hinweg (statt serieller
Ausführung wie bei BUSCO) sowohl möglich als auch sinnvoll.

## 2026-09-02 — Panel auf 5 Host-Repräsentanten reduziert (POC-Scope-Abweichung vom Dokument) + Repeat-Masking parallelisiert

**Wichtig:** Das Dokument (`multireferenzpanel_pav_workflow.md`) selbst
sieht KEINE Reduktion des 14-Genom-Katalogs vor — es geht durchgängig von
allen 14 Genomen aus. Die folgende Reduktion ist eine bewusste,
nutzergetriebene POC-Scope-Entscheidung (Rechenaufwand senken), keine
Vorgabe aus dem Dokument, und wird hier als Abweichung transparent
dokumentiert.

**Entscheidung:** Statt aller 14 Complete-Genome-Referenzen wird für
Repeat-Masking, Annotation und Panel-Bau nur noch **ein Repräsentant pro
Host-Typ** verwendet (5 Genome). `config/references.tsv` (aktiv, von
`workflow/Snakefile` gelesen) enthält jetzt nur diese 5 Zeilen; der volle
14-Genom-Katalog ist unverändert in
`config/references_full_catalog_14genomes.tsv` archiviert (die bereits
abgeschlossenen BUSCO-Ergebnisse für alle 14 bleiben gültig und werden
nicht verworfen, laufen aber nicht weiter durch nachfolgende
Pipeline-Schritte).

**Host-Typ-Tally und Korrektur einer Fehlklassifikation:** Ursprünglich
wurde `GCA036493215_1` als `host_group=unknown` geführt (das `host`-Feld
war im NCBI-BioSample-Datensatz leer). Eine gezielte Nachfrage bei der
NCBI Datasets API (`accession/GCA_036493215.1/dataset_report`) ergab
jedoch: Es handelt sich um **Br48**, ein Weizen-infizierendes Isolat aus
Brasilien (BioSample-Attribut `strain: "wheat infecting strain"`,
`geo_loc_name: Brazil`), UND die zugehörige BioProject-Beschreibung
lautet explizit "**Telomere-to-telomere genome assembly of Pyricularia
oryzae Br48**" (PRJDB14561) — ein selbstdeklariertes T2T-Assembly (7
Contigs = 7 Chromosomen, `contig_l50=3`). `host`/`host_group` wurden
entsprechend auf `Triticum aestivum`/`triticum` korrigiert. Damit sind es
nur **5 echte Host-Typen** unter den 14 Genomen (Oryza, Triticum,
Wildgrass/Lolium, Eleusine, Avena), nicht 6 wie zunächst angenommen.

**Repräsentantenauswahl:**

| Host-Typ | Kandidaten (14-Katalog) | Gewählt | Begründung |
|---|---|---|---|
| Oryza | 7015, Guy11, 95HPH4, 95085, P131 | **7015** | kanonischer Referenzstamm "70-15", Feldstandard in praktisch jeder vergleichenden *M.-oryzae*-Genomstudie |
| Triticum | GCA059330735_1, GCA059330115_1, GCA059330025_1, GCA059330365_1, Br48 | **GCA036493215_1 (Br48)** | einziges explizit T2T-deklariertes Genom im gesamten Panel — direkt relevant für saubere Starship-/Mini-Chromosom-Boundary-Calls (Abschnitt 6.3/8), wichtiger als der marginal höhere Contig-N50 der Alternativen |
| Wildgrass | LpKY97, GCA059329725_1 | **LpKY97** | etablierter, in der Literatur verwendeter Referenzstamm für die Lolium-Linie |
| Eleusine | GCA004346965_1 | GCA004346965_1 | einzige Complete-Genome-Option |
| Avena | GCA059329645_1 | GCA059329645_1 | einzige Complete-Genome-Option |

`GCA036493215_1` (Br48) behält seine Genome-ID (Accession-basiert statt
"Br48"), obwohl der Strain-Name jetzt bekannt ist — der bereits laufende
RepeatModeler-Pilotlauf (siehe oben) nutzt exakt diesen Wildcard-Wert;
eine Umbenennung hätte den bisherigen Fortschritt verwaist.

**Bekannter Trade-off:** Die Reduktion verliert Within-Host-Diversität
(z. B. 5 statt 1 Oryza-Genom, 5 statt 1 Triticum-Genom) — Starships, die
nur in einer Teilmenge der Isolate EINES Host-Typs vorkommen, werden vom
5-Genom-Panel nicht erfasst. Für die POC-Kernfrage (Nachweisbarkeit von
Starships/Accessory-Chromosomen ÜBER Host-Grenzen hinweg per Long-Read-
Mapping) ist das akzeptabel; für eine spätere Vollanalyse/Publikation
sollte auf den vollen 14-Genom-Katalog zurückgegriffen werden.

**Repeat-Masking parallelisiert:** Da RepeatModeler kaum RAM braucht
(siehe Pilotlauf-Befund oben), laufen jetzt alle 5 Panel-Genome
GLEICHZEITIG statt seriell: der bereits laufende Br48-Pilot (8 Threads,
unverändert weitergelaufen) plus ein zweiter, per `--nolock` parallel
gestarteter Snakemake-Lauf für die restlichen 4 Genome (`7015`,
`LpKY97`, `GCA004346965_1`, `GCA059329645_1`, je 3 Threads,
`config.yaml: threads_default` dafür von 8 auf 3 gesenkt). `--nolock`
ist hier sicher, da beide Läufe disjunkte Zielgenome und damit disjunkte
Ausgabedateien haben. Ressourcen-Check bei 5 parallelen Läufen: 4,5 GB
RAM (von 10 GB), Load Average ~10 (von 14 Kernen) — stabil, kein
OOM-Risiko.
