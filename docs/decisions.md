# Methodische Entscheidungen

## 2026-08-31 — Eindeutige Contig-IDs

**Entscheidung:** FASTA-Header werden zu
`<assembly_accession><separator><original_contig_id>` normalisiert.

**Begründung:** Starfish erwartet eine Genome-ID plus Feature-/Contig-ID.
Eindeutige IDs verhindern Kollisionen zwischen Isolaten.

**Konsequenz:** Zugehörige GFF-SeqIDs werden über dieselbe Mapping-Tabelle
synchron angepasst und anschließend gegen die FASTA validiert.

## 2026-08-31 — Trennzeichen `-` statt `_`

**Entscheidung:** Als Separator zwischen Isolat-ID und Contig-ID wird `-`
verwendet (`config/parameters.yaml: separator`), nicht `_` wie ursprünglich
im manuellen Starfish-Testlauf.

**Begründung:** NCBI-Assembly-Accessions wie `GCA_004346965.1` enthalten
selbst einen Unterstrich. Mit `_` als Separator ergab die zusammengesetzte
ID `GCA_004346965.1_CP034210.1` beim Aufsplitten `>2` statt `2` Komponenten;
Starfish konnte GenomeID und FeatureID nicht mehr eindeutig trennen
(Warnung: "is being parsed into >2 components using separator '_'").
Mit `-` (kommt weder in GCA-Accessions noch in RefSeq/GenBank-Contig-IDs
vor) verschwindet die Warnung vollständig, bei identischem YR-Ergebnis
(13 Kandidaten).

**Konsequenz:** `normalize_fasta_headers.py` bricht jetzt mit einer
Fehlermeldung ab, falls der gewählte Separator im Isolat-Namen selbst
vorkommt. `starfish annotate` wird konsequent mit `--separator` aufgerufen,
damit FASTA-Normalisierung und Starfish-Parsing konsistent bleiben.

## 2026-08-31 — GFF ist pro Isolat optional

**Entscheidung:** `config/samples.tsv` erlaubt eine leere `gff`-Spalte.
`normalize_gff_seqids` und `validate_fasta_gff_ids` laufen nur für Isolate,
die tatsächlich eine GFF-Datei angeben; `starfish_annotate_yr` läuft für
diese Isolate ohne `--gff` (rein de-novo via MetaEuk/HMM).

**Begründung:** Von 40 bisher heruntergeladenen *M. oryzae*-Assemblies hat
nur `GCA_004346965.1` (der Pilot) eine NCBI-Genannotation als GFF; alle
anderen liegen nur als FASTA + GenBank-Flatfile vor. Für die geplante
Skalierung auf viele Isolate ist das voraussichtlich der Regelfall, nicht
die Ausnahme. Ein Zwang zu vorhandener GFF hätte die Multi-Isolat-Pipeline
faktisch blockiert.

**Konsequenz:** Für den Drei-Isolat-Pilot wurden zwei zusätzliche,
wirtskontrastierende Isolate ohne GFF ergänzt: `GCA_004785725.2` (B71,
*Triticum aestivum*, Weizenblast-Referenzstamm) und `GCA_046718735.1`
(Guy11, *Oryza sativa*, meistgenutzter Laborstamm). Ergebnis: 13 (Pilot,
Fingerhirse), 12 (B71, Weizen) und 16 (Guy11, Reis) HMM-validierte
YR-Kandidaten, keine Header-Warnungen bei keinem der drei Isolate.

## 2026-08-31 — Mitochondriale Contigs von mChr-Kandidaten trennen

**Entscheidung:** `calculate_contig_metrics.py` bekommt eine `assembly_unit`-
Spalte, befüllt über das originale NCBI `sequence_report.jsonl` (Feld
`assemblyUnit`), verknüpft über die alte-ID/neue-ID-Mapping-Tabelle. Dafür
hat `config/samples.tsv` eine neue optionale Spalte `sequence_report`.

**Begründung:** In der ersten Contig-Metrik-Tabelle fiel bei zwei der drei
Pilotisolate je ein auffälliger Ausreißer-Contig auf: ca. 35 kb lang und mit
~28,6 % GC deutlich von allen anderen Contigs (1,9–8,8 Mb, ~48–51 % GC)
abweichend — in `GCA_004785725.2` (B71) `CP060338.1`, in `GCA_046718735.1`
(Guy11) `CM102890.1`. Prüfung im jeweiligen `sequence_report.jsonl` zeigt:
Beide sind als `"assemblyUnit": "non-nuclear"`,
`"assignedMoleculeLocationType": "Mitochondrion"` deklariert — es handelt
sich um das mitochondriale Genom, nicht um ein akzessorisches
Kern-Mini-Chromosom. Größe (~35 kb) und niedriger GC-Gehalt sind für
*M. oryzae*-mtDNA typisch. Der Pilot `GCA_004346965.1` enthält kein
separates MT-Contig in seiner Assembly.

**Konsequenz:** Jede spätere mChr-Klassifikation (`core_like`,
`accessory_candidate`, `mChr_candidate`, `uncertain`) muss Contigs mit
`assembly_unit == "non-nuclear"` vorab ausschließen. Ohne diesen Filter
wäre das Mitogenom fälschlich als starker mChr-Kandidat erschienen (klein,
Core-Gen-arm, stark abweichender GC-Gehalt) — ein Lehrbuchbeispiel für die
in `Starfish_und_MiniChromosomen_Analyseplan.md` (Abschnitt 11) genannte
Gefahr der Fehlzuordnung.

## 2026-08-31 — repeat_fraction über windowmasker, TE_fraction vertagt

**Entscheidung:** `repeat_fraction` wird über `windowmasker` (NCBI BLAST+,
bereits in `starfish_env` installiert) berechnet: Anteil weichmaskierter
(kleingeschriebener) Basen pro Contig nach dem Standard-Zweischritt-Ablauf
(`-mk_counts` → `-ustat ... -outfmt fasta`). Das ist eine
bibliotheksfreie, genomsequenz-basierte Repetitivitäts-Schätzung.

**Begründung:** Ein echtes `TE_fraction` (kuratierte Transposon-Familien,
z. B. LTR/DNA-Transposon-Klassen) erfordert RepeatMasker/EDTA mit einer
kuratierten, artspezifischen Repeat-Bibliothek — ein eigenständiges,
größeres Infrastrukturvorhaben. `windowmasker` liefert ohne zusätzliche
Bibliothek einen seriösen ersten Repetitivitäts-Proxy. Getestet: Pilot
zeigt 14,1 % genomweit maskierten Anteil, plausibel für *M. oryzae*.

**Konsequenz:** `TE_fraction` bleibt vorerst unbesetzt und wird erst
ergänzt, wenn eine RepeatMasker/EDTA-Umgebung mit passender Bibliothek
aufgesetzt wird. `repeat_fraction` darf nicht mit einer TE-Familienanalyse
verwechselt werden.

## 2026-08-31 — Erste mChr-Klassifikations-Heuristik

**Entscheidung:** `define_mchr_candidates.py` klassifiziert jeden
Kern-Contig (nicht `non-nuclear`) relativ zu einer pro Isolat berechneten
Core-Referenz (längengewichteter Mittelwert von GC und Repeat-Anteil sowie
Gendichte über Contigs ≥ `core_min_length_bp`, Default 1 Mb):

- `core_like`: Länge ≥ 1 Mb.
- Contigs ≤ 500 kb erhalten je 1 Evidenzpunkt für: Repeat-Anreicherung
  ≥ 0,10 über Referenz, GC-Abweichung ≥ 0,03 absolut, Gendichte ≤ 50 % der
  Referenz (nur wenn `gene_count` verfügbar ist).
  - ≥2 Evidenzlinien → `mChr_candidate`
  - genau 1 → `accessory_candidate`
  - 0 → `uncertain`
- Größe zwischen 500 kb und 1 Mb, oder keine Core-Referenz im Isolat
  auffindbar → `uncertain`.

Alle Schwellen liegen in `config/parameters.yaml` (`mchr:`-Sektion).

**Begründung:** Entspricht der in
`Starfish_und_MiniChromosomen_Analyseplan.md` (C2/C3) geforderten
Kombination mehrerer unabhängiger Evidenzlinien statt einer einzelnen
Schwelle, mit den aktuell tatsächlich verfügbaren Metriken.

**Wichtiger Befund/Caveat:** Im Drei-Isolat-Testlauf wurden bei
`GCA_046718735.1` (Guy11) fünf kleine Contigs (35–244 kb) als
`accessory_candidate`/`mChr_candidate` markiert (hohe `repeat_fraction`
von 33–96 %). Alle fünf sind laut NCBI `sequence_report.jsonl` als
`role: unplaced-scaffold` gekennzeichnet — Sequenzstücke, die der
Assembler wegen ihrer Repetitivität keinem der 7 Chromosomen zuordnen
konnte. Das ist eine ebenso plausible Erklärung wie ein echtes
akzessorisches Element (Fragmentierungs-/Fehlzuordnungsgefahr, siehe
Analyseplan Abschnitt 11). Deshalb wurde eine `role`-Spalte
(`assembled-molecule` vs. `unplaced-scaffold`, aus `sequence_report.jsonl`)
in die Ausgabetabelle aufgenommen, aber bewusst NICHT in die
Klassifikationslogik eingebaut — sie dient der manuellen Einordnung, nicht
einer automatischen Ausschlussregel. `GCA_004785725.2` (B71) zeigt keine
solchen Kandidaten (alle 8 Kern-Contigs `core_like`).

**Konsequenz:** Jede inhaltliche Interpretation eines
`mChr_candidate`/`accessory_candidate` muss die `role`-Spalte prüfen.
Ein `unplaced-scaffold`-Kandidat braucht zusätzliche, unabhängige Evidenz
(z. B. Synteny-Vergleich zu anderen Isolaten oder Read-Coverage), bevor er
ernsthaft als akzessorisches Chromosom in Betracht gezogen wird.

## 2026-08-31 — mChr-Größenschwellen an publizierte *M. oryzae*-Biologie angepasst

**Entscheidung:** `core_min_length_bp` wurde von 1 Mb auf 4 Mb angehoben,
`small_max_length_bp` von 500 kb auf 3 Mb. Contigs zwischen 3 und 4 Mb
sind bewusst eine Grauzone (`uncertain`), keine automatische Zuordnung.

**Begründung (Nutzerhinweis):** Mini-Chromosomen bei *M. oryzae* gelten als
supernumerär, typischerweise repeatreicher und genärmer als Core-Chromosomen;
publizierte Größen reichen von einigen hundert kb bis etwa 3 Mb. Die
vorherigen Schwellen (Core ab 1 Mb) hätten Contigs in genau diesem
publizierten mChr-Größenbereich automatisch als `core_like` durchgewunken,
ohne die übrigen Evidenzlinien überhaupt zu prüfen.

**Beobachtete Auswirkung im Drei-Isolat-Testlauf:** `GCA_004785725.2`
(B71) hat einen 1,9-Mb-Contig (`CP060337.1`), der vorher `core_like` war
und jetzt `uncertain` ist (GC- und Repeat-Abweichung vom Core-Referenzwert
liegen jeweils knapp unter der Nachweisschwelle — kein starkes Signal,
aber auch keine automatische Core-Zuordnung mehr). `GCA_046718735.1`
(Guy11) hat einen 3,9-Mb-Contig (`CM102889.1`), der jetzt in die 3–4-Mb-
Grauzone fällt. Beide verdienen manuelle Prüfung, bevor sie als Core oder
mChr gelten.

**Konsequenz:** Die Klassifikation ist konservativer geworden (mehr
`uncertain`-Fälle statt vorschneller `core_like`-Zuordnung). Das ist
gewünscht: ein unsicheres, dokumentiertes Ergebnis ist besser als eine
falsche Sicherheit (vgl. Analyseplan Abschnitt 14).

## 2026-08-31 — Telomer-Vollständigkeit und 5-Klassen-Schema

**Entscheidung:** `calculate_contig_metrics.py` sucht in den ersten/letzten
`telomere_window_bp` (Default 1000 bp) jedes Contigs nach einem
Tandem-Lauf von mindestens `telomere_min_repeats` (Default 5) Kopien des
kanonischen fungalen Telomer-Repeats `(TTAGGG)n` (oder seines
Revers-Komplements `(CCCTAA)n`, da die Ausrichtung des Contigs relativ zum
Chromosomenende vorab unbekannt ist). Ergebnis: `telomere_start`/
`telomere_end` (bool) pro Contig.

`define_mchr_candidates.py` nutzt das für eine 5. Klasse:
- `high_confidence_mChr`: ≥2 Evidenzlinien UND Telomere an beiden Enden.
- `mChr_candidate`: ≥2 Evidenzlinien, aber Telomere fehlen/unvollständig.
- `accessory_region` (vorher `accessory_candidate`, umbenannt auf
  Nutzer-Terminologie): genau 1 Evidenzlinie.
- `uncertain`, `core_like`, `excluded_non_nuclear`: unverändert.

**Begründung (Nutzerhinweis):** Long-Read-Assemblies liefern besonders
wertvolle mChr-Evidenz, wenn ein Contig vollständig von Telomer zu Telomer
reicht (Long Reads können repetitive Bereiche überspannen). Ein
kompletter, eigenständiger chromosomaler Contig mit Telomeren an beiden
Enden ist ein stärkeres Signal als Größe/GC/Repeat/Gendichte allein.

**Validierung an echten Daten:** Vor dem Einbau am Pilot-Isolat getestet:
echte Telomere zeigten 18–31 Tandem-Kopien, Nicht-Telomer-Enden 0–1 — klar
trennbar, `min_repeats=5` liegt sicher dazwischen. Im Drei-Isolat-Lauf
erreicht kein Contig `high_confidence_mChr` (keiner der Kandidaten hat
Telomere an beiden Enden); mehrere Core-Chromosomen selbst fehlt an einem
Ende das Telomer (z. B. `GCA_004346965.1-CP034204.1` an beiden, mehrere
weitere an einem) — plausibel für reale Assemblies (rDNA-Arrays u. Ä.
blockieren manche Chromosomenenden) und bestätigt, dass die Erkennung
nicht trivial immer "True" liefert.

**Noch nicht umgesetzt (aus dem vollständigen 5-Klassen-Schema):**
Fragmentierungsprüfung (ein zusammenhängender Contig vs. mehrere
unverbundene Stücke desselben mChr) und Core-Synteny-Vergleich (breites
kollineares Alignment zu Core-Chromosomen) fehlen noch — `core_like`
basiert aktuell nur auf Größe, nicht auf Syntenie. Diese sind der
nächste geplante Ausbauschritt und erfordern einen Aligner
(minimap2/nucmer), der noch nicht eingebunden ist.

## 2026-08-31 — Core-Synteny per minimap2-Selbst-Alignment

**Entscheidung:** Neue Regeln `align_genome_self` (minimap2, Preset
`asm5`, `--secondary=yes -N 5`, bereits in `starfish_env` vorhanden) und
`compute_core_synteny.py` berechnen pro Contig `core_synteny_coverage`:
den Anteil seiner Länge, der mit einem ANDEREN, core-großen Contig
desselben Isolats kollinear ist (gemergte Alignment-Intervalle, triviale
Self-Hits ausgeschlossen).

`define_mchr_candidates.py` nutzt das zweifach:
- **Veto:** `core_synteny_coverage ≥ 0.5` erzwingt `uncertain`, unabhängig
  von anderer Evidenz (breite Kollinearität zu einem Core-Chromosom
  spricht für Assembly-Fragment/Duplikat, nicht für ein eigenständiges
  Element).
- **Voraussetzung für `high_confidence_mChr`:** zusätzlich zu ≥2
  Evidenzlinien und Telomeren an beiden Enden muss `core_synteny_coverage
  < 0.1` sein (Kriterium "geringe Core-Syntenie" aus dem Nutzerschema).

**Begründung:** Ergänzt das letzte fehlende Kriterium aus dem
Nutzer-Schema für `high_confidence_mChr`. Ohne Synteny-Check hätte ein
großes, aber schlicht dupliziertes Stück eines Core-Chromosoms fälschlich
als eigenständiger mChr-Kandidat durchgehen können.

**Beobachtung im Drei-Isolat-Testlauf:** Bei `GCA_046718735.1` (Guy11)
richten sich die fünf `unplaced-scaffold`-Contigs (JAWCTR...) fast
ausschließlich **untereinander** aus, nicht gegen die 7 echten
Core-Chromosomen (`core_synteny_coverage = 0.0` für alle). Das bestätigt
"geringe Core-Syntenie" — deutet aber gleichzeitig an, dass diese fünf
Scaffolds möglicherweise redundante/überlappende Assemblierungsversuche
derselben repetitiven Region sind, nicht fünf unabhängige Elemente. Diese
Mutual-Redundanz zwischen NICHT-Core-Contigs wird von
`compute_core_synteny.py` bewusst nicht erkannt (nur Kollinearität zu
core-großen Contigs zählt) — das ist die separate, noch offene
Fragmentierungsprüfung aus dem Nutzerschema.

**Ergebnis:** Kein Contig im aktuellen Datensatz erreicht
`high_confidence_mChr` (keiner erfüllt gleichzeitig Evidenz + beidseitige
Telomere + geringe Syntenie). `GCA_046718735.1-JAWCTR010000013.1` bleibt
`mChr_candidate` (Syntenie niedrig bestätigt, aber keine Telomere).

**Noch offen:** Fragmentierungsprüfung (mehrere unverbundene, aber
zusammengehörige Contigs desselben mChr über gegenseitige Homologie
erkennen) ist die letzte fehlende Komponente aus dem 5-Klassen-Schema.

## 2026-08-31 — Multi-Referenzpanel für Core/Synteny-Konsens

**Entscheidung:** `config/references.yaml` definiert ein Referenzpanel mit
zwei Qualitätsstufen:
- `structural` (Chromosomen-nah/Long-Read, für Struktur/Synteny/mChr):
  70-15 (GCA_000002495.2, einzige mit kuratierter NCBI-Annotation), Guy11
  (GCA_046718735.1), B71 (GCA_004785725.2), CD156 (GCA_900474475.3,
  Eleusine-Linie), US71 (GCA_900474175.3, Setaria-Linie).
- `supplementary` (älter/fragmentiert, nur zur Diversitätserweiterung,
  nicht für strukturelle Platzierung): BdMeh16 (GCA_001675605.1,
  Bangladesch-MoT), EI9411 (GCA_001548775.1) und EI9604
  (GCA_001548785.1, beide zweite Eleusine-Sublinie).

"PH14 oder GY11" (dritte, global kontrastierende Reislinie) wurde auf
Nutzerwunsch ausgelassen, da keine Accession lokal auffindbar war und
keine geraten werden sollte.

`workflow/scripts/list_reference_core_contigs.py` bestimmt pro Referenz
die Core-Contigs (Länge ≥ `core_min_length_bp`, nicht-nukleär
ausgeschlossen, dokumentierte Akzessorik manuell ausgeschlossen via
`known_accessory_contigs`). Wichtig: `role` aus dem NCBI-Report ist bei
nicht-skaffoldierten Contig-Level-Referenzen (CD156, US71) für JEDEN
Contig `unplaced-scaffold`, unabhängig von der Größe — dort zählt nur die
Länge.

`workflow/scripts/compute_reference_synteny.py` richtet jedes Isolat
gegen jede strukturelle Referenz aus (minimap2, Preset `asm20`,
`--secondary=no`) und zählt `reference_core_hits` (Anzahl Referenzen mit
Coverage ≥ 0,3 gegen deren Core-Contigs). `reference_core_consensus`
(≥ `consensus_min_hits`, Default 2) überschreibt die interne
Größen-Heuristik zu `core_like` — eine mit mehreren unabhängigen
Referenzgenomen kollineare Sequenz gilt als Teil des konservierten
Core-Genoms, unabhängig vom Aussehen im eigenen Isolat. Bestätigte 0
Treffer zählen zusätzlich als Evidenzlinie für Akzessorik/mChr.

**Zwei Bugs während der Entwicklung gefunden und behoben (wichtig für
Ergebnisinterpretation):**

1. **Selbst-Vergleich-Kontamination:** B71 und Guy11 sind sowohl
   Test-Isolate als auch Panel-Referenzen. Ein Isolat gegen sich selbst
   als Referenz zeigt trivial 100 % Identität und hätte den eigenen
   Konsens künstlich aufgebläht. Fix: `structural_reference_ids_for_sample()`
   im Snakefile schließt eine Referenz aus, deren `accession` mit der
   `assembly_accession` des Samples übereinstimmt.

2. **Fehlender Mapping-Qualitäts-/Ziel-Konzentrations-Filter:** B71s
   1,9-Mb-Contig (`CP060337.1`, vermutetes dokumentiertes
   Mini-Chromosom) zeigte gegen 70-15/CD156/Guy11 hunderte kurze,
   über viele verschiedene Zielchromosomen verstreute Alignments mit
   überwiegend niedriger `mapq` — das Muster eines im Core-Genom
   mehrfach vorkommenden Repeat-/Transposon-Elements, nicht echter
   1:1-Syntenie. Die ursprüngliche Coverage-Berechnung summierte das über
   alle Ziel-Contigs hinweg auf und erzeugte so `reference_core_hits=3`
   (fälschlich `core_like`). Fix: Coverage wird jetzt als **beste
   Einzelziel-Coverage** (nicht Summe über alle Ziele) berechnet, und
   Alignments mit `mapq < 30` werden ausgeschlossen. Nach dem Fix:
   `reference_core_hits=1`, `consensus=False` → `uncertain` (weder
   sicher Core noch sicher akzessorisch — angemessen vorsichtig).

**Konsequenz:** Alle künftigen Synteny-Vergleiche (Self- und
Referenz-basiert) müssen mapq-gefiltert und ziel-konzentriert bewertet
werden; reine Coverage-Summen über mehrere Ziel-Contigs sind für
Repeat-reiche Sequenzen nicht aussagekräftig. B71s mutmaßliches
Mini-Chromosom bleibt nach Korrektur `uncertain`, nicht bestätigt oder
widerlegt — die Diskrepanz zu "B71 hat ein dokumentiertes
Mini-Chromosom" ist noch nicht mit der Originalpublikation
abgeglichen.

## 2026-08-31 — Starfish über `conda run -n starfish_env`

**Entscheidung:** Die Snakemake-Regel `starfish_annotate_yr` ruft Starfish
über `conda run -n starfish_env starfish annotate ...` auf, statt über
Snakemakes `conda:`-Direktive mit eigener Environment-Datei.

**Begründung:** `starfish_env` existiert bereits manuell eingerichtet
(inkl. Referenzdatenbanken unter `$CONDA_PREFIX/db`) und wird nicht von
Snakemake verwaltet. Ein Wechsel auf eine Snakemake-verwaltete
Environment-Definition ist ein separater, größerer Schritt (Datenbanken
müssten mitverwaltet werden) und ist für den aktuellen Pilotlauf nicht
nötig.

**Konsequenz:** Die restlichen Regeln (`normalize_fasta_headers`,
`normalize_gff_seqids`, `validate_fasta_gff_ids`) laufen weiterhin über
Snakemakes reguläre `conda:`-Direktive mit `workflow/envs/python.yaml`.
