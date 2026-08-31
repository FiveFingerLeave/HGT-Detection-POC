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
