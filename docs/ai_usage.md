# KI-Nutzungsprotokoll

## Grundsatz
KI dient für Erklärung, Codeentwürfe, Dokumentation und Review.
Alle wissenschaftlichen Entscheidungen, Parameter, Tests und Interpretationen
werden von mir geprüft und dokumentiert.

## 2026-09-01
- Verwendung: Vollständiger Projekt-Reset auf Basis von
  `Dokumentation/POC Workflow – HGT Machbarkeitsstudie M. oryzae.md`.
  Entfernung des vorherigen, methodisch anderen Pipeline-Ansatzes
  (Assembly-interne Contig-Klassifikation). Neuaufbau der Repo-Struktur,
  vier Snakemake-Regeldateien, vier Python-Skripte, R-Power-Analyse-Skript
  gemäß Dokument.
- Verwendung: Abruf des vollständigen NCBI-SRA-Katalogs für Pyricularia
  oryzae (3.902 Läufe) via E-Utilities, gefiltert nach WGS/Plattform.
- Verifikation: Alle neuen Python-Skripte mit pytest getestet (13 Tests).
  `validate_clustering_synthetic.py` tatsächlich ausgeführt (nicht nur
  geschrieben) — Ergebnis unten dokumentiert.
- Von der KI gefundenes Problem, selbst verifiziert: Ein WSL-Absturz
  (OOM-Kill, dokumentiert in `docs/decisions.md`-Vorgängerversion in der
  Git-Historie) hatte zwei Dokumentationsdateien
  (`Starfish_bisherige_Schritte.md`,
  `Starfish_und_MiniChromosomen_Analyseplan.md`) von der Festplatte
  entfernt, ohne dass dies beabsichtigt war. Aus der Git-Historie
  wiederhergestellt (Commit `dd8adf5`), vor dem endgültigen Reset-Commit
  geprüft.
- Von der KI gefundenes Problem, selbst verifiziert: `.gitignore` enthielt
  eine pauschale `data/`-Ausschlussregel, die auch die neu erstellten
  NCBI-Isolat-/Sequenzlisten von der Versionierung ausgeschlossen hätte.
  Vor dem Commit korrigiert (gezielte Ausschlüsse für große Sequenzdateien
  statt des gesamten Ordners).
- Ergebnis mit wichtigem Caveat: Der synthetische Clustering-Validierungstest
  (POC-Dokument Abschnitt 3) zeigt bei realistischer Monte-Carlo-Mittelung
  (20 Wiederholungen), dass ein einzelnes injiziertes HGT-Signal unter
  wenigen Kandidatenregionen nur in ~30 % der Fälle einen klaren
  ARI-Abfall erzeugt — nicht robust. Dies wurde transparent dokumentiert
  statt durch günstige Parameterwahl verdeckt (siehe `docs/decisions.md`).
