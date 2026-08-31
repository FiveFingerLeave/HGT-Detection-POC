# KI-Nutzungsprotokoll

## Grundsatz
KI dient für Erklärung, Codeentwürfe, Dokumentation und Review.
Alle wissenschaftlichen Entscheidungen, Parameter, Tests und Interpretationen
werden von mir geprüft und dokumentiert.

## 2026-08-31
- Verwendung: Aufbau der Repository-Grundstruktur (Snakemake/Conda-Projekt),
  Entwurf von `normalize_fasta_headers.py`, `normalize_gff_seqids.py` und
  `validate_fasta_gff_ids.py` samt pytest-Tests, Verdrahtung als
  Snakemake-Regeln (`workflow/rules/headers.smk`, `workflow/rules/starfish.smk`).
- Verifikation: Alle Skripte per pytest getestet (9 Tests, lokal ausgeführt);
  komplette Pipeline zweimal end-to-end mit dem echten Pilot-Isolat
  `GCA_004346965.1` laufen lassen (`snakemake --cores 6`), Logs geprüft.
- Ergebnis: 13 HMM-validierte YR-Kandidaten in `GCA_004346965.1`
  (identisch zum vorherigen manuellen Testlauf).
- Interpretation: Als YR_candidate dokumentiert, nicht als bestätigtes
  Starship (siehe `docs/decisions.md`).
- Von der KI gefundenes Problem, selbst verifiziert: Der ursprünglich
  verwendete Separator `_` kollidierte mit dem Unterstrich in
  NCBI-Accessions (`GCA_004346965.1`), was Starfish-Warnungen
  ("is being parsed into >2 components") auslöste. Vor der Änderung im
  Log bestätigt, nach der Änderung (`separator: "-"`) verifiziert, dass
  die Warnungen verschwinden und die YR-Kandidatenzahl unverändert bleibt.
