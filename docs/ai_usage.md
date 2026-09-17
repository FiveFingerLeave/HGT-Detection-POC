# AI Usage Log

## Principle
AI is used for explanation, code drafts, documentation, and review.
All scientific decisions, parameters, tests, and interpretations
are reviewed and documented by me.

## 2026-09-01
- Use: Complete project reset based on
  `Documentation/POC_Workflow_HGT_Feasibility_Study_M_oryzae.md`.
  Removal of the previous, methodologically different pipeline approach
  (assembly-internal contig classification). Rebuilt the repo structure,
  four Snakemake rule files, four Python scripts, R power-analysis script
  per the document.
- Use: Retrieval of the complete NCBI SRA catalog for Pyricularia
  oryzae (3,902 runs) via E-Utilities, filtered by WGS/platform.
- Verification: All new Python scripts tested with pytest (13 tests).
  `validate_clustering_synthetic.py` actually executed (not just
  written) — result documented below.
- Problem found by the AI, self-verified: A WSL crash
  (OOM kill, documented in an earlier version of `docs/decisions.md` in the
  git history) had removed two documentation files
  (`Starfish_Progress_Log.md`,
  `Starfish_and_MiniChromosome_Analysis_Plan.md`) from disk,
  unintentionally. Restored from the git history (commit `dd8adf5`),
  verified before the final reset commit.
- Problem found by the AI, self-verified: `.gitignore` contained
  a blanket `data/` exclusion rule that would also have excluded the
  newly created NCBI isolate/sequence lists from version control.
  Corrected before the commit (targeted exclusions for large sequence files
  instead of the whole folder).
- Result with an important caveat: The synthetic clustering validation test
  (POC document section 3) shows, under realistic Monte Carlo averaging
  (20 repetitions), that a single injected HGT signal among few
  candidate regions produces a clear ARI drop in only ~30 % of cases —
  not robust. This was documented transparently instead of being
  concealed through a favorable choice of parameters (see `docs/decisions.md`).
