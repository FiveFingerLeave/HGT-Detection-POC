# Barragan – Magnaporthe oryzae Starship/HGT-Pipeline

Reproduzierbare Snakemake-Pipeline für Starfish-, Starship-, Mini-Chromosomen- und Isolatvergleichsanalysen bei *Magnaporthe oryzae*.

## Voraussetzungen

- Conda/Mamba
- Environment `barragan-workflow` (Snakemake, Python 3.12, pandas, biopython, pytest, ruff):

```bash
mamba env create -f workflow/envs/python.yaml
conda activate barragan-workflow
```

## Ordnerübersicht

- `config/` – Sample-Liste, Parameter und Referenzpfade
- `workflow/` – Snakefile, Regeln (`rules/`), Skripte (`scripts/`), Conda-Environments (`envs/`)
- `input/` – Manifeste und Metadaten (kein Git für Rohdaten)
- `results/`, `logs/`, `benchmarks/` – generierte Outputs (nicht versioniert)
- `docs/` – Methoden, Entscheidungen, KI-Nutzungsprotokoll
- `tests/` – pytest-Tests für eigene Skripte

## Testlauf (Dry-run)

```bash
snakemake --snakefile workflow/Snakefile --cores 6 --use-conda --dry-run --printshellcmds
```

Details zum Setup: [Dokumentation/Entwicklungsumgebung_und_KI_Integration.md](Dokumentation/Entwicklungsumgebung_und_KI_Integration.md)
