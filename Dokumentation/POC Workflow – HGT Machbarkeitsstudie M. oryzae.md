# POC Workflow: Machbarkeitsstudie HGT in *M. oryzae*

Guideline für die technische Umsetzung in VS Code · Bezug: Dissertationsexposé „Population Genomics of Horizontal Gene Transfer in Magnaporthe oryzae", WP1.2–WP1.4, Q1–Q3 · Stand 2026\-09\-01

## Zweck dieses Dokuments

Dieses Dokument ist eine **Arbeitsgrundlage für VS Code**, kein Endergebnis. Es definiert einen kleinskaligen Proof\-of\-Concept (POC), der vor Antragstellung/Förderzusage zeigen soll, dass die im Exposé beschriebene Methodik (WP1.2–WP1.4) technisch funktioniert und die geplante Auswertung (Q3, GLMM) statistisch sinnvoll ist. Es deckt die vier angefragten Punkte ab:

1. Funktionalität des Multireferenzpanels
2. Anzahl Starships und weiterer HGT\-Kandidaten im POC
3. Clustering von Core\- vs. Kandidatenregionen nach Host Type
4. Power\-Analyse für das geplante GLMM (Q3)

**Wichtige Annahme:** In dieser Session liegen keine realen Sequenzdaten vor (nur das Exposé\-PDF im Projekt). Die Schritte 1–3 sind daher als **ausführbares Protokoll** mit konkreten Tools/Befehlen/Parametern formuliert, bereit zur Anwendung, sobald Referenzgenome und ein Isolat\-Subset (z. B. aus dem 417\-Isolat\-Datensatz, Barragan et al. 2024, Suppl. Table S5) verfügbar sind. Schritt 4 (Power\-Analyse) ist datenunabhängig und kann bereits jetzt mit angenommenen Effektgrößen simuliert werden — diese Annahmen sind explizit markiert und müssen aktualisiert werden, sobald erste Pilotdaten aus WP1.3 vorliegen.

## Empfohlene Repo\-Struktur

```text
poc-hgt-magnaporthe/
├── data/
│   ├── references/            # Multireferenzpanel (FASTA + GFF), s. Abschnitt 1
│   └── isolates_poc/          # Illumina/Nanopore-Subset für den POC
├── envs/                      # conda/mamba environment.yml je Tool
├── workflow/
│   ├── Snakefile              # oder nextflow main.nf
│   ├── rules/                 # mapping.smk, pav.smk, starfish.smk, cluster.smk
│   └── scripts/
├── results/
│   ├── mapping/
│   ├── pav_calls/
│   ├── starfish/
│   └── clustering/
├── power_analysis/
│   └── glmm_power_sim.R       # s. Abschnitt 4
└── README.md                  # ← dieses Dokument (als Export)
```

Ein Workflow\-Manager (Snakemake oder Nextflow) wird empfohlen, da WP1.3 explizit eine „modulare Pipeline" mit kalibrierbaren Schwellenwerten fordert — reproduzierbar und später direkt auf den vollen 417\-Isolat\-Datensatz sowie das HPC\-Cluster der CAU skalierbar.

* * *

## 1\. Funktionalität Multireferenzpanel

### Ziel

Prüfen, ob das im Exposé beschriebene Mapping\-Verfahren (Kurzreads gegen ein Multireferenzpanel, coverage\-/PAV\-basierte Kandidatendetektion) technisch funktioniert und die vorgeschlagenen Schwellenwerte (80–90 % Breadth bei ≥5× Tiefe \= „present", \<10–20 % Breadth \= „absent") kalibrierbar sind.

### Referenzauswahl (POC\-Umfang: 4–6 Genome)

Kriterium laut Exposé: Abdeckung mehrerer Host\-Lineages plus möglichst hochwertige, publizierte Long\-Read\-Assemblies.

| Lineage / Host | Referenz | Quelle |
| --- | --- | --- |
| Reis (*Oryza sativa*) | 70\-15 (near\-complete) | Cheng et al. 2025 |
| Weizen (*Triticum aestivum*) | Wheat\-blast\-Referenzisolat | Yoshida et al. 2016 |
| Finger Hirse (*Eleusine coracana*) | Eleusine\-Lineage\-Isolat | Chiapello et al. 2015 |
| Wildgras\-Lineages | ≥1 Isolat aus Barragan et al. 2022/2024 | CAU / bioRxiv\-Supplement |
| *M. grisea* (Schwestertaxon, Negativkontrolle) | Referenzassembly | Chiapello et al. 2015 |

### Pipeline\-Schritte

```bash
# 1. QC der Kurzreads
fastp -i isolate_R1.fq.gz -I isolate_R2.fq.gz \
      -o qc_R1.fq.gz -O qc_R2.fq.gz --json qc.json

# 2. Panel indizieren (ein Index pro Referenz, oder kombiniertes Panel)
bwa-mem2 index references/panel_combined.fa
# Alternative: minimap2 -d panel.mmi references/panel_combined.fa

# 3. Mapping
bwa-mem2 mem -t 8 references/panel_combined.fa qc_R1.fq.gz qc_R2.fq.gz \
  | samtools sort -@4 -o mapping/isolate.sorted.bam
samtools index mapping/isolate.sorted.bam

# 4. Coverage / Breadth pro Kandidatenregion
mosdepth --by candidate_regions.bed -t 4 pav_calls/isolate mapping/isolate.sorted.bam

# 5. PAV-Call (Skript, Schwellenwerte parametrisiert)
python scripts/pav_call.py \
  --regions candidate_regions.bed \
  --mosdepth pav_calls/isolate.regions.bed.gz \
  --present-breadth 0.85 --present-depth 5 \
  --absent-breadth 0.15 \
  --out pav_calls/isolate.pav.tsv
```

### Kalibrierung & Erfolgskriterien

- Kalibrierung der Schwellenwerte an **bekannten Mini\-Chromosomen\-Regionen** (mChrA, Barragan et al. 2024) und Replikatproben.
- Ziel: **Sensitivität ≥ 90 %, Spezifität ≥ 95 %**, finale Schwellenwerte über F1\-Score\-Maximierung (Grid\-Search über Breadth\-/Tiefe\-Kombinationen).
- Sensitivitätsanalyse: Schwellenwerte in 5\-%\-Schritten variieren, PAV\-Calls gegen long\-read\-bestätigte Regionen validieren.

### POC\-Machbarkeitsfragen (zu beantworten)

- [ ] Mappen Reads aus divergenten Lineages (z. B. Wildgras\-Isolat gegen Reis\-Referenz) noch zuverlässig genug für ein Coverage\-Signal?
- [ ] Wie stark beeinflusst Referenzbias (welche Referenz im Panel „gewinnt") die PAV\-Calls? → Test: gleiches Isolat gegen Einzelreferenz vs. kombiniertes Panel.
- [ ] Skaliert der Workflow rechnerisch auf den vollen Datensatz (Laufzeit/Speicher pro Isolat × 417)?

* * *

## 2\. Starship\- und HGT\-Kandidatenzählung im POC

### Ziel

Für den POC\-Subset (s. o., 4–6 Referenzen \+ 10–20 Isolate) ermitteln, **wie viele Starships und weitere HGT\-Kandidatenregionen** mit den beiden im Exposé beschriebenen komplementären Methoden gefunden werden.

### Zwei komplementäre Detektionswege (laut Exposé, Methode a\+b)

**a) Coverage\-basiert (Multireferenzpanel, s. Abschnitt 1):** PAV\-Regionen aus `pav_call.py`, gefiltert auf Kandidaten mit diskordantem Clustering (s. Abschnitt 3).

**b) Struktur\-basiert (Starfish / Stargraph):**

```bash
# Starfish: Detektion von Starship-Elementen anhand Captain-Tyrosin-Rekombinase
# + charakteristischer Flanking-Signaturen (Gluck-Thaler & Vogan 2024)
starfish annotate -T 8 -x isolate_id \
  -a references/isolate_assembly.fa \
  -g references/isolate_assembly.gff3 \
  -o results/starfish/isolate_id

starfish insert -x isolate_id \
  -a results/starfish/isolate_id.starships.bed \
  -o results/starfish/isolate_id.inserts

# Stargraph: Pan-genomischer Graph-Ansatz zur Konsistenzprüfung über Isolate hinweg
stargraph build --assemblies references/*.fa --starships results/starfish/*.bed \
  -o results/starfish/pangenome_graph
```

### Kreuzvalidierung (laut Exposé)

- Coverage\-basierte PAV\-Kandidaten → auf Starship\-Hallmarks screenen (Captain\-Protein, terminale Repeats, Target\-Site\-Duplikationen).
- Starfish\-Kandidaten → auf Coverage\-Diskordanz über Lineages prüfen.
- Klassifikation jeder Kandidatenregion in **Vektorklasse**\: `mini-chromosome`, `Starship`, `beides`, `unklar`.

### Erwarteter POC\-Output

Eine Kandidatentabelle als zentrales Zwischenergebnis:

| Feld | Beschreibung |
| --- | --- |
| `region_id` | Eindeutige ID der Kandidatenregion |
| `vector_class` | mChr / Starship / beides / unklar |
| `size_kb` | Größe (Starships: 25–700 kb laut Definition) |
| `n_isolates_present` | Anzahl POC\-Isolate mit „present"\-Call |
| `host_lineages_present` | Betroffene Host\-Lineages (für Abschnitt 3) |
| `captain_detected` | ja/nein (nur Starships) |
| `evidence_tier` | hochkonfident / explorativ (laut Risikobewertung Q1) |

**Zielgröße für den POC:** Da echte Zahlen erst nach Lauf der Pipeline vorliegen, ist die POC\-Erfolgsmetrik nicht „N Starships", sondern: *Wird mindestens 1 bekannte, publizierte Starship\-Region (O'Donnell et al. 2025) korrekt wiedergefunden, und liefert die Pipeline eine nicht\-triviale, aber handhabbare Kandidatenzahl (grobe Erwartung im niedrigen zweistelligen Bereich für \~10–20 Isolate, keine feste Vorgabe)?* Dies als Sensitivitäts\-Sanity\-Check, nicht als Preprint\-Zahl.

* * *

## 3\. Clustering Core\- vs. Kandidatenregionen nach Host Type

### Ziel

Die im Exposé gezeigte Logik aus **Abbildung 2** reproduzieren: Core\-Genom\-Marker sollen Isolate nach Host\-assoziierter Lineage clustern, Kandidaten\-HGT\-Regionen sollen dagegen **diskordantes** (lineage\-übergreifendes) Clustering zeigen — das zentrale Evidenzkriterium für HGT (Q1).

### Methodik

1. **Core\-Genom\-Matrix:** PAV\- oder SNP\-Matrix aus core\-genomweiten Markern (außerhalb der Kandidatenregionen).
2. **Kandidaten\-Matrix:** PAV\-Matrix nur aus den in Abschnitt 2 identifizierten Kandidatenregionen.
3. Für beide Matrizen: PCA (oder MDS), gefärbt nach `host_lineage` (Reis, Weizen, Eleusine, Wildgras).
4. **Diskordanz quantifizieren:** Cluster aus PCA (k\-means oder hierarchisches Clustering, k \= Anzahl bekannter Lineages) mit den bekannten Host\-Labels vergleichen via **Adjusted Rand Index (ARI)**.
   - Erwartung Core\-Matrix: **hoher ARI** (Cluster ≈ Host\-Lineage).
   - Erwartung Kandidaten\-Matrix: **niedrigerer ARI**, mit klar identifizierbaren lineage\-übergreifenden Gruppen (\= HGT\-Signatur).

```python
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

def cluster_and_score(pav_matrix: pd.DataFrame, host_labels: pd.Series, k: int):
    pca = PCA(n_components=min(10, pav_matrix.shape[1]))
    coords = pca.fit_transform(pav_matrix.values)
    clusters = KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(coords)
    ari = adjusted_rand_score(host_labels, clusters)
    return coords, clusters, ari

core_coords, core_clusters, core_ari = cluster_and_score(core_pav, host_labels, k=4)
cand_coords, cand_clusters, cand_ari = cluster_and_score(candidate_pav, host_labels, k=4)
print(f"Core ARI: {core_ari:.2f}  |  Kandidaten-ARI: {cand_ari:.2f}")
```

### Empfehlung: Methodenvalidierung an synthetischen Daten (jetzt schon möglich)

Bevor echte Daten vorliegen, kann die Analyselogik selbst an einem synthetischen Toy\-Datensatz validiert werden (bekannte „Ground Truth" HGT\-Isolate künstlich eingebaut, Erwartung: Pipeline muss injizierte Diskordanz korrekt erkennen). Ein Skript dafür (`scripts/validate_clustering_synthetic.py`) ist als eigenständiger Test empfohlen, der prüft:

- Core\-Matrix ohne injizierte HGT\-Events → ARI nahe 1.
- Kandidaten\-Matrix mit künstlich eingefügtem cross\-lineage sharing (z. B. 15 % der Isolate teilen eine „fremde" Kandidatenregion) → ARI signifikant niedriger, betroffene Isolate im Kandidaten\-PCA sichtbar gruppiert.

Dies dient als **Unit\-Test der Analysemethode**, nicht als biologisches Ergebnis.

### Interpretationskriterien

- [ ] Core\-ARI deutlich \> Kandidaten\-ARI (qualitativ, kein fixer Schwellenwert im Exposé vorgegeben)
- [ ] Visuelle Diskordanz reproduziert Muster aus Abb. 2 (Wildgras\-/Reis\-Isolate clustern in Kandidatenregionen zusammen)
- [ ] Ergebnis konsistent über mehrere k\-Werte / Clustering\-Methoden (Robustheitscheck)

* * *

## 4\. Power\-Analyse für das GLMM (Q3)

### Ziel

Abschätzen, ob eine POC\-/Pilot\-Stichprobe (bzw. die spätere volle Feldstudie) ausreichend Power hat, um die im Exposé spezifizierten Effekte zu detektieren: Poisson\- bzw. Negativ\-Binomial\-GLMM mit fixen Effekten (Donor\-Präsenz, Recipient\-Präsenz, Vektortyp, Donor×Vektor\- und Recipient×Vektor\-Interaktion), Random Intercepts für Feld und Jahr, Offset \= log(Anzahl beprobter Isolate pro Feld\-Jahr).

### Empfohlenes Vorgehen: simulationsbasierte Power\-Analyse (Monte\-Carlo)

**Empfohlenes Tooling:** R mit `glmmTMB` (bzw. `lme4::glmer.nb`) \+ `simr` — Standardwerkzeuge für genau dieses Modell in der Ökologie/Evolutionsbiologie (vgl. Bolker et al. 2009, im Exposé zitiert). *Hinweis: In dieser Cloud\-Session ist kein R installiert — dieser Abschnitt liefert ein direkt in VS Code lauffähiges Codegerüst; als Fallback ohne R eignet sich `statsmodels.genmod.bayes_mixed_glm.PoissonBayesMixedGLM` in Python, mit geringerer Genauigkeit bei kleinen Clusterzahlen.*

### POC\-Stichprobenszenarien (Annahmen — mit Pilotdaten aus WP1.3 zu aktualisieren)

| Szenario | Felder | Jahre | Isolate / Feld\-Jahr | Gesamt\-Isolate | Angenommene HGT\-Baseline\-Rate | Angenommener Effekt Donor×Vektor (Rate Ratio) |
| --- | --- | --- | --- | --- | --- | --- |
| Konservativ (Minimal\-POC) | 5 | 2 | 10 | 100 | 5 % | 1\.5× |
| Moderat (realistischer POC) | 10 | 2 | 20 | 400 | 8 % | 1\.8× |
| Optimistisch (erweiterter POC) | 15 | 3 | 30 | 1\.350 | 12 % | 2\.0× |
| Volle Studie (Referenz) | \>20 | ≥3 | variabel | \~1.000\+ | unbekannt | — |

Diese Werte sind **Platzhalter**, keine publizierten Effektgrößen — HGT\-Ereignisraten sind laut Exposé selbst (Risikobewertung Q3) noch unbekannt und potenziell selten/ungleich verteilt.

### Codegerüst (`power_analysis/glmm_power_sim.R`)

```r
library(glmmTMB)
library(simr)

set.seed(1)

simulate_scenario <- function(n_fields, n_years, n_iso_per_fy,
                               baseline_rate, effect_rr) {
  df <- expand.grid(field = factor(1:n_fields), year = factor(1:n_years))
  df$log_n_isolates <- log(n_iso_per_fy)
  df$donor     <- rbinom(nrow(df), 1, 0.5)
  df$recipient <- rbinom(nrow(df), 1, 0.5)
  df$vector    <- factor(sample(c("mChr", "Starship", "both"),
                                 nrow(df), replace = TRUE))

  beta <- c(intercept = log(baseline_rate),
            donor = 0.3, recipient = 0.2,
            vectorStarship = 0.2, vectorboth = 0.4,
            donor_vectorStarship = log(effect_rr))

  # linearer Prädiktor inkl. Random Intercepts (Feld, Jahr) + Poisson-Ziehung
  field_re <- rnorm(n_fields, 0, 0.4)[df$field]
  year_re  <- rnorm(n_years, 0, 0.3)[df$year]
  lp <- beta["intercept"] + beta["donor"] * df$donor +
        beta["recipient"] * df$recipient +
        ifelse(df$vector == "Starship", beta["vectorStarship"],
               ifelse(df$vector == "both", beta["vectorboth"], 0)) +
        ifelse(df$donor == 1 & df$vector == "Starship",
               beta["donor_vectorStarship"], 0) +
        field_re + year_re + df$log_n_isolates

  df$hgt_count <- rpois(nrow(df), exp(lp))
  df
}

power_for_scenario <- function(..., nsim = 200) {
  df <- simulate_scenario(...)
  fit <- glmmTMB(hgt_count ~ donor * vector + recipient * vector +
                   offset(log_n_isolates) + (1 | field) + (1 | year),
                 family = poisson, data = df)
  powerSim(fit, test = fixed("donor:vectorStarship"), nsim = nsim)
}

# Szenarien aus obiger Tabelle durchlaufen und Power (%) protokollieren
scenarios <- list(
  konservativ  = list(n_fields = 5,  n_years = 2, n_iso_per_fy = 10, baseline_rate = 0.05, effect_rr = 1.5),
  moderat      = list(n_fields = 10, n_years = 2, n_iso_per_fy = 20, baseline_rate = 0.08, effect_rr = 1.8),
  optimistisch = list(n_fields = 15, n_years = 3, n_iso_per_fy = 30, baseline_rate = 0.12, effect_rr = 2.0)
)

results <- lapply(scenarios, function(s) do.call(power_for_scenario, s))
```

### Interpretation & Machbarkeitsentscheidung

- **Faustregel:** Bei so wenigen Clustern (5–15 Felder × 2–3 Jahre) ist die Power für **Interaktionseffekte** (Donor×Vektor) typischerweise deutlich geringer als für Haupteffekte — dies ist der kritischste Test in der Simulation.
- Wenn selbst das „optimistische" POC\-Szenario \< 80 % Power für den Interaktionseffekt zeigt: (a) Modell vereinfachen (Interaktion nur explorativ, wie im Exposé als Fallback vorgesehen: „stepwise, beginning with simple models"), oder (b) Stichprobenumfang für die volle Studie entsprechend dimensionieren.
- Sobald WP1.3 erste reale HGT\-Kandidatenraten liefert, `baseline_rate` und `effect_rr` durch beobachtete Werte ersetzen und Simulation wiederholen — die jetzige Analyse ist ein **Startpunkt**, kein finaler Beleg.

### Machbarkeits\-Checkliste

- [ ] Power ≥ 80 % für Haupteffekte (Donor, Recipient, Vektortyp) im moderaten Szenario?
- [ ] Power ≥ 80 % für mindestens einen Interaktionseffekt in mind. einem realistischen Szenario?
- [ ] Sensitivität der Power gegenüber Random\-Effect\-Varianz (Feld/Jahr) geprüft?
- [ ] Vergleich Poisson vs. Negativ\-Binomial (Überdispersion) eingeplant, falls reale Zähldaten überdispers sind?

* * *

## Zusammenfassung: Offene Punkte für die Umsetzung in VS Code

- [ ] Zugriff auf Referenzgenome (Abschnitt 1) und Isolat\-Subset organisieren (HPC\-Cluster CAU)
- [ ] Environments aufsetzen: `bwa-mem2`/`minimap2`, `samtools`, `mosdepth`, `starfish`, `stargraph`, Python (`pandas`, `scikit-learn`), R (`glmmTMB`, `simr`)
- [ ] Snakemake/Nextflow\-Grundgerüst gemäß Repo\-Struktur oben anlegen
- [ ] Synthetischen Clustering\-Validierungstest (Abschnitt 3) implementieren und laufen lassen — unabhängig von echten Daten sofort möglich
- [ ] Power\-Simulation (Abschnitt 4) in R lokal/auf HPC ausführen und Szenariotabelle mit Ergebnissen füllen
- [ ] Nach erstem Pipeline\-Lauf: Kandidatentabelle (Abschnitt 2) und ARI\-Werte (Abschnitt 3) in dieses Dokument nachtragen, Annahmen in Abschnitt 4 aktualisieren
