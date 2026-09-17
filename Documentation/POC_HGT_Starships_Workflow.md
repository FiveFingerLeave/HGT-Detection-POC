# Proof of Concept: Reliable Detection of HGT Candidates and Starships

## 1. Objective (adjusted)

**Not:** Exact determination of HGT frequency.
**Rather:** Demonstrating that the workflow reliably identifies candidates for horizontal gene transfer (HGT) and Starships — plus determining a viable sample size and model complexity for the later GLMM-based HGT analysis in the dissertation.

The POC delivers three types of results:
1. **Methodological validation** — the reference panel and the mapping/detection workflow reliably and reproducibly identify accessory regions/Starships.
2. **Descriptive level** — frequency and variability of Starships/accessory regions across clonal lineages and host plants.
3. **HGT evidence** — concrete candidate cases (nearly identical Starships in distantly related clonal lineages, divergent insertion sites, phylogenetic incongruence) as a basis for later long-read verification.

Additional secondary goal: simulation-based recommendation for genome number and GLMM complexity in the main study.

---

## 2. Phase 1 — Reference Panel Construction

- Selection of 20–30 isolates from as many different clonal lineages and host plants as possible.
- Priority: already available long-read assemblies (NCBI search for assembly level "Complete Genome"/"Chromosome", sequencing technology PacBio/Nanopore).
- For gaps without public long-read data: own long-read sequencing, otherwise exclusion from the core panel.
- **QC:** BUSCO (completeness), assembly statistics (N50, contig count), contamination check (e.g. BlobTools).
- **Annotation:** consistent gene prediction (e.g. funannotate/BRAKER) — important, since later ortholog comparisons rely on consistent gene models.
- **Repeat masking:** RepeatModeler + RepeatMasker (Starships are located in repetitive/TE-rich regions, hence essential for clean boundary calls).

**Output:** curated panel of N genomes with consistent annotation + QC report.

---

## 3. Phase 2 — Starship Identification in the Reference Panel

- Detection with an established Starship finder (e.g. *starfish*: Captain gene search, boundary calling, cargo extraction).
- Manual curation of the calls: boundary plausibility, target site duplications, Captain family assignment.
- Construction of a Starship catalog per genome: ID, coordinates, Captain family, cargo genes, insertion site/locus, flanking sequence.

**Output:** Starship catalog (table) + sequence/feature files per element.

---

## 4. Phase 3 — Presence/Absence Matrix (Pan-Accessory Genome)

- Pangenome tool (Panaroo/PIRATE/Roary) for gene PAV across all panel genomes.
- Additional synteny-/alignment-based detection of non-genic accessory regions (e.g. via whole-genome alignment/pangenome graph).
- Linking the PAV matrix with the Starship catalog: which accessory blocks are Starships, which originate from other sources (e.g. other TEs, deletions)?

**Output:** binary PAV matrix (isolates × elements), annotated as Starship vs. non-Starship.

---

## 5. Phase 4 — Saturation/Rarefaction Analysis (Core of the POC)

- Random, repeated subsampling of increasing genome numbers (1 … N) from the panel, many permutations (100–1000 repetitions per size) — analogous to gene accumulation curves in the pangenome field.
- Tools: R (`vegan::specaccum`, `micropan`) or a custom script.
- Plot cumulative number of newly discovered Starships/accessory regions against genome number.
- Fit a saturation model (Heaps' law or asymptotic model) → check whether and at which genome number the curve flattens (e.g. < 5% new elements per additional genome).
- The result provides a direct justification: "panel size X is sufficient" or "further reference genomes needed".

**Output:** saturation curve + model parameters + panel size recommendation.

---

## 6. Phase 5 — Validation on Test Isolates (Leave-one-out + Real Test Isolates)

- **Leave-one-out validation:** remove each reference genome from the panel once, re-call PAV using a "short-read-like" mapping approach (reads/simulated reads against the remaining panel), compare with known long-read ground truth → sensitivity, specificity, precision of the method.
- **Real test isolates:** map additional isolates (short-read data only) against the reference panel, PAV calling via coverage in Starship boundary regions, compare results with the panel catalog.

**Output:** validation metrics (sensitivity/specificity) as evidence of method reliability.

---

## 7. Phase 6 — Phylogenetic Incongruence (HGT Evidence)

- Core genome phylogeny from single-copy orthologs (IQ-TREE/RAxML) as a reference tree (= clonal lineage relatedness).
- Separate phylogeny of the Starship Captain genes/cargo genes across the panel.
- Topology comparison (Robinson-Foulds distance, ALE/Notung reconciliation) → strong incongruence = HGT candidate.
- Additional evidence: nearly identical Starships (>99% identity) in phylogenetically distant clonal lineages; different insertion sites of the same/similar Starship in different lineages.

**Output:** list of prioritized HGT candidates with supporting evidence, as a basis for later long-read confirmation.

---

## 8. Phase 7 — Simulation: Sample Size & Model Complexity for the Later GLMM

- Simulation of synthetic PAV datasets under controlled conditions (known HGT rate, number of clonal lineages, genome number, effect size) — e.g. via simulated phylogenies (R: `ape`, `phytools`) with simulated transfer events.
- Fitting GLMMs of increasing complexity (fixed effects only → + random effect clonal lineage → + random effect host plant → interaction terms) at varying sample sizes (e.g. N = 10/20/30/50).
- Evaluation: convergence rate, bias/variance of the estimators, power to detect effects, stability across bootstrap repetitions.
- Result: recommendation for the minimum genome number and maximum sensible model complexity for the main study.

**Output:** simulation report with power/bias curves and a concrete sample size/model recommendation.

---

## 9. Phase 8 — Synthesis & Go/No-Go

Consolidation into a POC report with four key statements:
1. Is the panel size sufficient? (from the saturation analysis)
2. How variable/frequent are Starships/accessory regions? (descriptive)
3. Which concrete HGT candidates were found, with what strength of evidence?
4. Which sample size/model complexity is recommended for the dissertation?

Go/No-Go criteria, e.g.: saturation curve flattens AND validation sensitivity > defined threshold (e.g. 90%) → the workflow is viable, scaling to the full cohort is justified.

---

## 10. Tool Stack (Overview)

| Step | Tools |
|---|---|
| Assembly/QC | hifiasm/Flye, BUSCO, BlobTools |
| Annotation | funannotate / BRAKER |
| Repeat masking | RepeatModeler, RepeatMasker |
| Starship detection | starfish |
| Pangenome/PAV | Panaroo, PIRATE, Roary |
| Genome alignment | minimap2, Cactus/pggb |
| Rarefaction | R (vegan, micropan) |
| Phylogeny | IQ-TREE/RAxML, ALE/Notung |
| Simulation & GLMM | R (ape, phytools, lme4/glmmTMB) |

---

## 11. Pipeline Structure for an Executing Agent

```
poc_hgt_starships/
├── 00_data/           # Raw genomes, NCBI downloads, metadata (clonal lineage, host)
├── 01_assembly_qc/     # Assembly, BUSCO, contamination check
├── 02_annotation/       # Gene models, repeat masking
├── 03_starship_calls/    # starfish output, curated catalog
├── 04_pav_matrix/       # Pangenome PAV, Starship annotation of the matrix
├── 05_saturation/       # Subsampling scripts, rarefaction curves, model fit
├── 06_validation/       # Leave-one-out, test isolate mapping, metrics
├── 07_phylogeny/        # Core tree, Starship gene trees, incongruence tests
├── 08_simulation_glmm/   # Simulation scripts, GLMM fits, power analysis
└── 09_report/          # Summary POC report
```

**Execution order:** 01 → 02 → 03 → 04 → (05 in parallel with 06) → 07 → 08 → 09.
Each folder should contain its own README with input/output specification and a reproducible script (Snakemake or Nextflow rule), so that the agent can run the entire workflow in a version-controlled and repeatable manner.
