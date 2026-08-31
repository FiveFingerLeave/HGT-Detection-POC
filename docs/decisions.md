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
