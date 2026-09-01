# Proof of Concept: Zuverlässige Detektion von HGT-Kandidaten und Starships

## 1. Zielsetzung (angepasst)

**Nicht:** Exakte Bestimmung der HGT-Häufigkeit.
**Sondern:** Nachweis, dass der Workflow zuverlässig Kandidaten für horizontalen Gentransfer (HGT) und Starships identifiziert — plus Bestimmung einer tragfähigen Stichprobengröße und Modellkomplexität für die spätere GLMM-basierte HGT-Analyse in der Dissertation.

Das POC liefert drei Ergebnistypen:
1. **Methodische Validierung** — das Referenzpanel und der Mapping-/Detektions-Workflow erkennen akzessorische Regionen/Starships zuverlässig und reproduzierbar.
2. **Deskriptive Ebene** — Häufigkeit und Variabilität von Starships/akzessorischen Regionen über Klonlinien und Wirtspflanzen hinweg.
3. **HGT-Indizien** — konkrete Kandidatenfälle (nahezu identische Starships in entfernten Klonlinien, abweichende Insertionsstellen, Phylogenie-Inkongruenz) als Grundlage für spätere Long-Read-Verifikation.

Zusätzliches Nebenziel: Simulationsbasierte Empfehlung für Genomzahl und GLMM-Komplexität in der Hauptstudie.

---

## 2. Phase 1 — Referenzpanel-Aufbau

- Auswahl von 20–30 Isolaten aus möglichst vielen unterschiedlichen Klonlinien und Wirtspflanzen.
- Priorität: bereits verfügbare Long-Read-Assemblies (NCBI-Suche nach Assembly-Level „Complete Genome"/„Chromosome", Sequencing-Technologie PacBio/Nanopore).
- Für Lücken ohne öffentliche Long-Read-Daten: eigene Long-Read-Sequenzierung, sonst Ausschluss aus dem Kernpanel.
- **QC:** BUSCO (Vollständigkeit), Assembly-Statistiken (N50, Kontiganzahl), Kontaminationscheck (z. B. BlobTools).
- **Annotation:** einheitliche Genvorhersage (z. B. funannotate/BRAKER) — wichtig, da spätere Ortholog-Vergleiche auf konsistenten Genmodellen beruhen.
- **Repeat-Masking:** RepeatModeler + RepeatMasker (Starships liegen in repetitiven/TE-reichen Regionen, daher essenziell für saubere Boundary-Calls).

**Output:** kuratiertes Panel von N Genomen mit einheitlicher Annotation + QC-Report.

---

## 3. Phase 2 — Starship-Identifikation im Referenzpanel

- Detektion mit einem etablierten Starship-Finder (z. B. *starfish*: Captain-Gen-Suche, Boundary-Calling, Cargo-Extraktion).
- Manuelle Kuration der Calls: Boundary-Plausibilität, Target-Site-Duplications, Captain-Familienzuordnung.
- Aufbau eines Starship-Katalogs pro Genom: ID, Koordinaten, Captain-Familie, Cargo-Gene, Insertionsort/-locus, flankierende Sequenz.

**Output:** Starship-Katalog (Tabelle) + Sequenz-/Feature-Dateien pro Element.

---

## 4. Phase 3 — Presence/Absence-Matrix (Pan-Accessory-Genom)

- Pangenom-Tool (Panaroo/PIRATE/Roary) für Gen-PAV über alle Panel-Genome.
- Ergänzend synteny-/alignment-basierte Detektion nicht-genischer akzessorischer Regionen (z. B. via Whole-Genome-Alignment/Pangenom-Graph).
- Verknüpfung der PAV-Matrix mit dem Starship-Katalog: welche akzessorischen Blöcke sind Starships, welche stammen aus anderen Quellen (z. B. andere TEs, Deletionen)?

**Output:** binäre PAV-Matrix (Isolate × Elemente), annotiert nach Starship vs. Nicht-Starship.

---

## 5. Phase 4 — Sättigungs-/Rarefaktionsanalyse (Kern des POC)

- Zufälliges, wiederholtes Subsampling steigender Genomzahlen (1 … N) aus dem Panel, viele Permutationen (100–1000 Wiederholungen je Größe) — analog zu Gen-Akkumulationskurven im Pangenom-Bereich.
- Werkzeuge: R (`vegan::specaccum`, `micropan`) oder eigenes Skript.
- Kumulative Zahl neu entdeckter Starships/akzessorischer Regionen gegen Genomzahl plotten.
- Sättigungsmodell fitten (Heaps'-Law oder asymptotisches Modell) → Prüfen, ob und ab welcher Genomzahl die Kurve abflacht (z. B. < 5 % neue Elemente pro zusätzlichem Genom).
- Ergebnis liefert direkte Begründung: „Panelgröße X ist ausreichend" oder „weitere Referenzgenome nötig".

**Output:** Sättigungskurve + Modellparameter + Empfehlung zur Panelgröße.

---

## 6. Phase 5 — Validierung an Testisolaten (Leave-one-out + echte Testisolate)

- **Leave-one-out-Validierung:** jedes Referenzgenom einmal aus dem Panel entfernen, mit „Short-Read-artigem" Mapping-Ansatz (Reads/simulierte Reads gegen Restpanel) PAV neu callen, mit bekanntem Long-Read-Ground-Truth vergleichen → Sensitivität, Spezifität, Präzision der Methode.
- **Echte Testisolate:** zusätzliche Isolate (nur Short-Read-Daten) gegen das Referenzpanel mappen, PAV-Calling über Coverage in Starship-Boundary-Regionen, Ergebnisse mit dem Panel-Katalog abgleichen.

**Output:** Validierungsmetriken (Sensitivität/Spezifität) als Nachweis der Methodenzuverlässigkeit.

---

## 7. Phase 6 — Phylogenetische Inkongruenz (HGT-Indizien)

- Core-Genom-Phylogenie aus Single-Copy-Orthologen (IQ-TREE/RAxML) als Referenzbaum (= Klonlinien-Verwandtschaft).
- Separate Phylogenie der Starship-Captain-Gene/Cargo-Gene über das Panel.
- Topologievergleich (Robinson-Foulds-Distanz, ALE/Notung-Rekonziliation) → starke Inkongruenz = HGT-Kandidat.
- Zusätzliche Indizien: nahezu identische Starships (>99 % Identität) in phylogenetisch entfernten Klonlinien; unterschiedliche Insertionsstellen desselben/ähnlichen Starships in verschiedenen Linien.

**Output:** Liste priorisierter HGT-Kandidaten mit Belegen, als Grundlage für spätere Long-Read-Bestätigung.

---

## 8. Phase 7 — Simulation: Stichprobengröße & Modellkomplexität für spätere GLMM

- Simulation synthetischer PAV-Datensätze unter kontrollierten Bedingungen (bekannte HGT-Rate, Anzahl Klonlinien, Genomzahl, Effektgröße) — z. B. über simulierte Phylogenien (R: `ape`, `phytools`) mit simulierten Transferereignissen.
- GLMMs steigender Komplexität fitten (nur Fixed Effects → + Random Effect Klonlinie → + Random Effect Wirtspflanze → Interaktionsterme) bei variierender Stichprobengröße (z. B. N = 10/20/30/50).
- Bewertung: Konvergenzrate, Bias/Varianz der Schätzer, Power zur Effektdetektion, Stabilität über Bootstrap-Wiederholungen.
- Ergebnis: Empfehlung für minimale Genomzahl und maximal sinnvolle Modellkomplexität für die Hauptstudie.

**Output:** Simulationsreport mit Power-/Bias-Kurven und konkreter Stichproben-/Modellempfehlung.

---

## 9. Phase 8 — Synthese & Go/No-Go

Zusammenführung zu einem POC-Bericht mit vier Kernaussagen:
1. Panelgröße ausreichend? (aus Sättigungsanalyse)
2. Wie variabel/häufig sind Starships/akzessorische Regionen? (deskriptiv)
3. Welche konkreten HGT-Kandidaten wurden gefunden, mit welcher Evidenzstärke?
4. Welche Stichprobengröße/Modellkomplexität wird für die Dissertation empfohlen?

Go/No-Go-Kriterien z. B.: Sättigungskurve flacht ab UND Validierungssensitivität > definierter Schwelle (z. B. 90 %) → Workflow ist tragfähig, Skalierung auf volle Kohorte gerechtfertigt.

---

## 10. Werkzeug-Stack (Übersicht)

| Schritt | Tools |
|---|---|
| Assembly/QC | hifiasm/Flye, BUSCO, BlobTools |
| Annotation | funannotate / BRAKER |
| Repeat-Masking | RepeatModeler, RepeatMasker |
| Starship-Detektion | starfish |
| Pangenom/PAV | Panaroo, PIRATE, Roary |
| Genom-Alignment | minimap2, Cactus/pggb |
| Rarefaktion | R (vegan, micropan) |
| Phylogenie | IQ-TREE/RAxML, ALE/Notung |
| Simulation & GLMM | R (ape, phytools, lme4/glmmTMB) |

---

## 11. Pipeline-Struktur für einen ausführenden Agenten

```
poc_hgt_starships/
├── 00_data/           # Rohgenome, NCBI-Downloads, Metadaten (Klonlinie, Wirt)
├── 01_assembly_qc/     # Assembly, BUSCO, Kontaminationscheck
├── 02_annotation/       # Genmodelle, Repeat-Masking
├── 03_starship_calls/    # starfish-Output, kuratierter Katalog
├── 04_pav_matrix/       # Pangenom-PAV, Starship-Annotation der Matrix
├── 05_saturation/       # Subsampling-Skripte, Rarefaktionskurven, Modellfit
├── 06_validation/       # Leave-one-out, Testisolat-Mapping, Metriken
├── 07_phylogeny/        # Core-Baum, Starship-Gen-Bäume, Inkongruenz-Tests
├── 08_simulation_glmm/   # Simulationsskripte, GLMM-Fits, Power-Analyse
└── 09_report/          # Zusammenfassender POC-Bericht
```

**Ausführungsreihenfolge:** 01 → 02 → 03 → 04 → (05 parallel zu 06) → 07 → 08 → 09.
Jeder Ordner sollte ein eigenes README mit Input/Output-Spezifikation und ein reproduzierbares Skript (Snakemake- oder Nextflow-Regel) enthalten, damit der Agent den gesamten Workflow versionskontrolliert und wiederholbar ausführen kann.
