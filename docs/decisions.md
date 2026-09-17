# Methodological Decisions

## 2026-09-01 — Project reset based on the POC workflow document

**Decision:** The entire previous pipeline setup (assembly-based contig
classification: `normalize_fasta_headers.py`, Starfish integration,
repeat/telomere/synteny heuristics, multi-reference panel for core
consensus, BUSCO/QUAST QC) was removed. Rebuilt from scratch based on
`Documentation/POC_Workflow_HGT_Feasibility_Study_M_oryzae.md`, which is
tied to the dissertation proposal "Population Genomics of Horizontal Gene
Transfer in Magnaporthe oryzae" (WP1.2–WP1.4).

**Rationale:** The POC document describes a fundamentally different
methodological approach from what had previously been built up
iteratively: short-read mapping against a multi-reference panel with
coverage-/PAV-based candidate detection (bwa-mem2/mosdepth), complemented
by Starfish/Stargraph (structure-based) and PCA/ARI clustering
(host-lineage discordance as HGT evidence) — not the previously built
assembly-internal contig heuristic. Maintaining both approaches in
parallel in the same repo would have caused confusion (user request:
"only information relevant to the POC should be kept in the repo").

**What was preserved:**
- Git history (all previous commits remain retrievable)
- `Documentation/` (all planning notes, including the new POC document)
- `assemblies/` (4.1 GB of already downloaded NCBI genomes, not tracked
  in Git, a purely local cache)
- NCBI isolate/sequence lists in `data/` and `input/metadata/`

**What was rebuilt from scratch (see below for details):**
- Repo structure exactly following the layout recommended in the POC
  document (`data/references/`, `data/isolates_poc/`, `envs/`,
  `workflow/`, `results/`, `power_analysis/`)
- Four environment definitions (`envs/mapping.yaml`, `starfish.yaml`,
  `python.yaml`, `r.yaml`)
- Snakemake skeleton (`workflow/Snakefile` + `rules/{mapping,pav,
  starfish,cluster}.smk`), closely modeled on the commands specified in
  the document, ready to run once `config/samples.tsv` and
  `data/references/panel_manifest.tsv` are populated
- Four new scripts: `pav_call.py`, `build_pav_matrix.py`,
  `cluster_and_score.py`, `validate_clustering_synthetic.py`
- `power_analysis/glmm_power_sim.R` (adopted from the document)

## 2026-09-01 — Comprehensive NCBI sequence dataset (SRA, not just assemblies)

**Decision:** In addition to the already existing assembly list
(`data/ncbi_m_oryzae_assemblies.tsv`/`.jsonl`, ~46 entries), the complete
NCBI SRA catalog for *Pyricularia oryzae* was retrieved
(`data/ncbi_m_oryzae_sra_runs.tsv`, via NCBI E-Utilities
esearch/efetch, `rettype=runinfo`): **3,902 sequencing runs**.
Filtered subsets:
- `data/ncbi_m_oryzae_sra_wgs_illumina.tsv` — 1,754 WGS short-read runs
  (Illumina) — candidate pool for the POC isolate subset (Section 1)
- `data/ncbi_m_oryzae_sra_wgs_longread.tsv` — 193 WGS long-read runs
  (Nanopore/PacBio) — candidate pool for reference-panel candidates

**Rationale:** The POC approach (short-read mapping against a
multi-reference panel) requires actual raw reads per isolate, not just
finished assemblies. The previous assembly list (46 entries) covers only
a fraction of the *M. oryzae* sequencing data actually available at
NCBI.

**Distribution by library strategy:** 1,971 WGS, 1,204 RNA-Seq, 435
ChIP-Seq, 86 WCS, 55 other, 43 WGA, 30 ncRNA-Seq, 25 bisulfite-seq, 20
RIP-Seq, 7 amplicon. By platform: 3,639 Illumina, 164 Oxford Nanopore,
51 PacBio, 19 DNBSEQ, 17 LS454, 12 BGISEQ.

**Consequence:** Before selecting the actual POC isolate subset (4–6
references + 10–20 isolates), isolates need to be selected from
`ncbi_m_oryzae_sra_wgs_illumina.tsv` according to host-lineage coverage
(rice, wheat, finger millet, wild grass depending on availability) —
not yet done; `config/samples.tsv` is currently only an empty template.

## 2026-09-01 — PAV calling: breadth from mosdepth thresholds, not just mean depth

**Decision:** `pav_call.py` computes breadth (the fraction of a region
with coverage ≥ `present_depth`) from mosdepth's `--thresholds` output
(`<prefix>.thresholds.bed.gz`), not from the mean depth in
`<prefix>.regions.bed.gz` alone.

**Rationale:** The POC document names both thresholds (breadth AND
depth) as criteria ("80–90% breadth at ≥5× depth"), but its example
invocation passes `pav_call.py` only the `regions.bed.gz` file (pure
mean depth per region) — no genuine breadth can be derived from that (a
region can have high mean depth from a small, extremely well-covered
sub-portion even though most of the region is not covered at all). The
`mosdepth --by` rule was therefore extended with
`--thresholds {present_depth}`.

**Consequence:** The `mosdepth_coverage` rule in `workflow/rules/pav.smk`
now produces two output files (`.regions.bed.gz` for depth,
`.thresholds.bed.gz` for breadth); `pav_call.py` accepts both separately
(`--mosdepth-regions`, `--mosdepth-thresholds`).

## 2026-09-01 — Region classification (core/candidate) via a shared BED file

**Decision:** `data/references/candidate_regions.bed` gets a 5th column,
`region_class` (`core` or `candidate`), instead of keeping core markers
and candidate regions in separate files/separate PAV runs.
`pav_call.py` passes this column through, and `build_pav_matrix.py`
then splits the combined `candidate_table.tsv` into two matrices for
Section 3 (PCA/ARI clustering).

**Rationale:** The document conceptually separates "core genome markers"
(genome-wide, outside the candidate regions) from "candidate regions"
(from Section 2), but does not specify a concrete file structure for
this. A shared BED file with a class column avoids a duplicated
mapping/mosdepth run (once per region class) and keeps the assignment
traceable in a single place.

## 2026-09-01 — Synthetic clustering validation test: an honest finding on method sensitivity

**Decision:** `validate_clustering_synthetic.py` (POC document Section 3,
"already possible now") runs as a Monte Carlo simulation over 20
repetitions by default, not as a single run.

**Rationale/finding:** A single test run with an injected cross-lineage
sharing of ONE candidate region (as literally described in the document:
"15% of isolates share a 'foreign' candidate region") shows strongly
different results depending on the random seed: with 5 candidate regions
and 30% recipient isolates, only ~30% of individual runs (seeds) showed
a clear ARI drop (>0.2); the averaged ARI drop over 20 repetitions is
only ~0.15 — not reliable enough for a robust diagnostic criterion. With
more candidate regions (e.g., 10–30), the signal disappears even more,
since unchanged, strictly lineage-faithful regions dominate the
PCA/k-means structure.

**Important finding for the actual methodology (not just a test
artifact):** An aggregated ARI value over the ENTIRE candidate matrix is
only weakly sensitive to a single genuine HGT signal when it is diluted
among several inconspicuous candidate regions. For the real analysis, a
**per-region discordance test** is therefore recommended in addition to
the aggregated ARI (e.g., lineage purity per individual candidate
region), before relying solely on an overall PCA/ARI metric. This is not
yet implemented.

**Consequence:** The feasibility question from Section 3 of the document
("the pipeline must correctly detect injected discordance") is **not
robustly met** by the current method as literally described in the
document — this needs to be communicated transparently before a funding
decision/application, not concealed through favorable parameter choices.

## 2026-09-01 — Per-region ARI implemented and directly compared against the aggregated method

**Decision:** `cluster_and_score.py` gets a new function
`per_region_ari()`: instead of a PCA over the entire candidate matrix,
the ARI between its presence/absence pattern and the known host lineage
is computed **for each candidate region individually** (missing/
"uncertain" values are excluded per region, not imputed).
`validate_clustering_synthetic.py` now runs both methods (A = aggregated
PCA/k-means ARI as described in the document, B = per-region ARI) on the
same synthetic data and compares them directly. The new output
`results/clustering/candidate_per_region_ari.tsv` is also available in
the real pipeline run (`workflow/rules/cluster.smk`).

**Result of the direct comparison (Monte Carlo, 20 repetitions):**

| Scenario | Affected isolates | Method A (aggregated) | Method B (per-region) |
|---|---|---|---|
| 5 candidate regions, 30% recipients | 12/40 | 30% detection rate | **95%** detection rate |
| 30 candidate regions, 15% recipients (literal document parameters) | 6/40 | **0%** detection rate | **70%** detection rate |

**Important secondary methodological finding:** With 4 host lineages but
only a binary (0/1) trait per region, the per-region ARI structurally
**cannot** be close to 1.0 — by the pigeonhole principle, at least two of
the four lineages must randomly share the same bit value, even with no
HGT at all (in the runs above: "clean" regions average ARI ≈ 0.38–0.39).
What is informative is therefore not the absolute value but the
**relative** distance: is the actually manipulated region clearly below
the average of the remaining regions? This is the case (0.17–0.26 vs.
0.38–0.39) and is measured via the detection rate ("correctly identified
as the most discordant region"), not via an absolute threshold.

**Consequence:** The per-region test is markedly more sensitive than the
aggregated approach in both tested scenarios, but even with the
document's own parameters (30 regions, only 6 affected isolates) reaches
only a 70% detection rate — also no guaranteed win. For the real POC run,
it is recommended to use `candidate_per_region_ari.tsv` as the primary
diagnostic tool and `ari_summary.tsv` (aggregated) only as a supplement,
not the other way around.

## 2026-09-01 — Per-region ARI operationalized: classification vs. plain ranking

**Decision:** `classify_region_discordance()` in `cluster_and_score.py`
turns the raw per-region ARI score into an actual per-region decision
(`discordant` / `concordant` / `uncertain`), via a robust modified
z-score (median + MAD, not mean + SD — so that individual genuine
outliers do not themselves distort the "normal" baseline). Output:
`results/clustering/candidate_per_region_ari.tsv` gets a
`classification` column; the threshold is configurable via
`config/parameters.yaml: clustering.discordance_z_threshold`.

**Calibration (Monte Carlo, standard scenario: 5 candidate regions, 30%
recipients, 20–30 repetitions):**

| z-threshold | Sensitivity | False-positive rate |
|---|---|---|
| 0.3 | 56.7% | 2.5% |
| 0.5 | 56.7–70.0% | 2.5% |
| 0.7 | 40.0% | 0.0% |
| 1.0 (original default) | 36.7–45.0% | 0.0% |

Default set to **z = 0.5** (best sensitivity/false-positive ratio
found).

**Important finding:** Even with a calibrated threshold, the
classification (at ~57–70% sensitivity) remains markedly behind the
plain ranking method ("is this region the lowest-scoring of all?", 95%
detection rate in the same scenario). Reason: with only ~5 candidate
regions per batch, the median/MAD estimate is itself statistically
unstable (small sample) — an absolute threshold cannot be calibrated as
reliably with so few regions as a simple rank comparison.

**Consequence for the actual POC run:** Given the realistically small
candidate counts expected for a POC (low double digits, per the
document's expectation), **sorting by per-region ARI (ascending)** is
the more reliable diagnostic tool, not the fixed classification
threshold. The `classification` column remains useful as a fast, coarse
filter with a very low false-positive rate, but should not serve as the
sole basis for decisions.

## 2026-09-01 — Three bug fixes during the first end-to-end run on real data

**Bug 1 (path):** `workflow/rules/regions.smk` had
`conda: "../envs/python.yaml"` (only one `../`), even though the rule
files live under `workflow/rules/` and the envs are in the project root
(`envs/python.yaml`) — what's needed is `"../../envs/python.yaml"`, as
in all other `.smk` files. As a result, the run aborted at
`core_marker_regions` with `WorkflowError: Error recording metadata`,
after all 6 Starfish runs had already completed. Fixed; the other three
`.smk` files already had the correct path.

**Bug 2 (resources, OOM):** bwa-mem2 was killed by the Linux OOM killer
while mapping `TH3` (`bwa-mem2.avx2` alone: 4.5 GB RSS). The WSL2 VM's
default limit (50% of host RAM = 7.5 GB of 16 GB) was too tight for
this — the same basic type of problem as the earlier BUSCO OOM. Since
the host actually has 16 GB of RAM, `C:\Users\<user>\.wslconfig` was
newly created (`memory=11GB`, `swap=4GB`) and WSL restarted
(`wsl --shutdown`). In addition, `map_to_panel`, `map_longread_to_panel`,
and `fastp_qc` in `mapping.smk` now get an explicit `threads:` value
(previously only `params.threads`, i.e., Snakemake reserved only 1 core
per job and, at `--cores 6`, could have allowed several 8-thread
bwa-mem2 processes to run simultaneously — with `threads:` they are now
correctly serialized instead).

**Bug 3 (missing NaN handling):** `cluster_and_score()` passed the PAV
matrix to `sklearn.decomposition.PCA` unchecked, which categorically
rejects NaN values ("uncertain" calls). This never surfaced in the
previous tests/synthetic validations, since those used only clean 0/1
data — but real coverage data routinely contains "uncertain" calls
(ambiguous breadth/depth). Fix: regions that are NaN for EVERY isolate
sample are dropped before the PCA (no information); remaining NaNs are
imputed with the column mean (presence frequency across the other
isolates) — a neutral default treatment for missing genotype-like calls.
A test was added
(`test_cluster_and_score_imputes_nan_instead_of_crashing`).

## 2026-09-01 — First end-to-end POC run on real data: results and limitations

**Data basis:** 6 reference genomes (70-15/Guy11 = rice, B71 = wheat,
CD156 = finger millet, US71 = foxtail millet, M_grisea = finger-millet
grass as a sister-taxon negative control), 5 short-read isolates
(MZ5-1-6, E34 = finger millet; TH3 = rice; U167 = finger-millet grass;
Py13.1.023 = wheat) + 1 long-read isolate (ZM1-2 = wheat, Nanopore). 92
Starship/YR candidate regions (Starfish `annotate`, all 6 references
successfully annotated, 12–18 YR loci per genome) + 120 core-marker
windows (20 per reference) in `data/references/candidate_regions.bed`.

**PAV calling result (`results/pav_calls/candidate_table.tsv`):**

| Region class | absent | present | uncertain |
|---|---|---|---|
| candidate (Starship) | 415 | 55 | 82 |
| core | 435 | 74 | 211 |

**Clustering (`results/clustering/ari_summary.tsv`):** Aggregated ARI
identical for both the core and the candidate matrix at 0.423 (k=4, n=6
isolates) — with only 6 isolates and 4 classes, the number of possible
partitions is so small that a coincidence carries little meaning; **this
metric is not reliable at the current POC sample size** (see also the
previously documented general weakness of aggregated ARI).

**Per-region ARI (`results/clustering/candidate_per_region_ari.tsv`):** 7
of 30 evaluable Starship regions classified as "discordant" (z=0.5):
`starship_B71_10`, `starship_CD156_3`, `starship_CD156_6`,
`starship_CD156_7`, `starship_CD156_12`, `starship_Guy11_3`,
`starship_M_grisea_10`. With only 6 isolates, this test is also
statistically very unstable (see the earlier calibration results at
n=40) — suitable as a shortlist for targeted manual follow-up, not as
evidence.

**Cross-lineage presence (manual supplementary analysis, not part of the
pipeline):** 13 cases in which a Starship candidate region from a
reference genome of a DIFFERENT host lineage than the isolate's own was
marked "present" (e.g., `E34`, a finger-millet isolate, shows presence
for Starship loci from `B71`/wheat, `US71`/foxtail millet, and
`M_grisea`/finger-millet grass) — this would be the expected picture of a
cross-lineage mobile element, but given the sample size it should be
treated as a hypothesis, not as proof. For the core windows there is
even more cross-lineage presence, with 34 cases, concentrated mainly in
`Py13.1.023` (21 of 34) — this suggests that the current core-marker
windows (randomly distributed windows on the largest contig, not
verified single-copy orthologs) capture genuine, biologically expected
conservation between references, are not reliably lineage-diagnostic,
and are additionally influenced by differences in sequencing depth
between isolates.

**Most important data limitation:** `U167` (finger-millet grass) has NO
"present" call in any of the 212 regions (neither core nor candidate) —
all 92 candidate and 120 core calls are "absent" or "uncertain", even for
regions from its own reference lineage (`M_grisea`). The raw data for
`U167` (SRR14705972), at 34 + 37 MB of compressed FASTQ, are noticeably
small compared to the other isolates (250 MB – 1.2 GB) — the sequencing
depth for this isolate is presumably fundamentally insufficient to
produce reliable presence calls against the 6-genome panel. `MZ5-1-6`
shows a similar, if less extreme, pattern (0 "present" calls for the
candidate regions). **A reliable statement requires either more
sequencing depth for these isolates or their exclusion from the
PAV-based analysis.**

**Overall conclusion:** The pipeline now runs completely end to end (QC →
mapping → PAV → clustering) and yields real, interpretable intermediate
results — but at n=6 isolates and uneven sequencing depth, neither the
aggregated ARI nor the per-region discordance list should be treated as
standalone evidence for HGT. This is expected for the feasibility study
(POC purpose: demonstrate pipeline functionality, not significance
already); for a funding application it should be explicitly stated that
sample size and depth differences currently limit statistical power.

## 2026-09-01 — R environment set up, power analysis run — result: 0% power in all three scenarios

**Background:** `power_analysis/glmm_power_sim.R` was carried over
verbatim from the POC document during the project reset, but never
executed (R was not installed). `envs/r.yaml` was now set up via
`conda env create -f envs/r.yaml` (r-base 4.3.3, glmmTMB 1.1.9, simr
1.0.7, lme4 1.1-37) and the script was run for the first time.

**Bug 1 (print output):** `summary(powerSim(...))` returns a
`data.frame` with columns `successes/trials/mean/lower/upper` — **not**
`Power`. The original snippet accessed `summary(x)$Power`, which always
returned `NULL`; `sprintf()` with a `NULL` argument returns
`character(0)` in R, causing the entire output loop to silently print
NOTHING (no error, no line) — the first run therefore falsely looked
like a clean but empty success. Fixed: now prints `s$mean` (× 100 for
percent), along with `s$successes`/`s$trials`.

**Bug 2 (methodological, deeper):** In the original snippet, `powerSim()`
is called without explicitly anchoring the effect to be tested
(`effect_rr`) in the model — it simply uses the coefficient estimated
from ONE randomly simulated dataset as the test target. `simr` itself
explicitly flags this as an **"observed power" calculation**
(`observedPowerWarning`) — a known biased, high-variance procedure, not
a genuine simulation-based power test for a SPECIFIED effect. The
actually correct `simr` workflow (`fixef(fit) <- ...`, to explicitly set
the target effect before simulation) was tested and fails on a genuine
compatibility gap: `simr` 1.0.7's `fixef<-` method dispatches to an S4
generic that does not recognize glmmTMB objects (`Error in
getClass(cl): "glmmTMB" is not a defined class`). A clean fix would
presumably require switching to `lme4::glmer`/`glmer.nb` (named as an
alternative in the POC document itself) — this was not implemented in
this session (effort/scope). Only the print bug was fixed; the result
below is the "observed power" approximation, not a rigorous power
analysis.

**Result (200 simulations per scenario):**

| Scenario | Fields×years×isolates/field-year | Power (donor×Starship interaction) | Successful sims |
|---|---|---|---|
| Conservative | 5×2×10 (100 isolates) | **0.0% (0.0–1.8%)** | 0/200 |
| Moderate | 10×2×20 (400 isolates) | **0.0% (0.0–1.8%)** | 0/200 |
| Optimistic | 15×3×30 (1,350 isolates) | **0.0% (0.0–1.8%)** | 0/200 |

In addition: in the "conservative" scenario, the interaction term
`donor:vectorStarship` itself — i.e., exactly the effect to be tested —
is dropped from the model due to rank deficiency (`dropping columns from
rank-deficient conditional model`). The design (only 5×2=10 field-year
cells) is too sparse to even estimate the full fixed-effects structure
(`donor*vector + recipient*vector`, 9 coefficients). In the
moderate/optimistic scenario the target term does remain in the model,
but power is nevertheless exactly 0% across all 200 repetitions (0 of
200 simulations yielded a significant result) — well below what would be
expected from random variation alone.

**Interpretation:** Because of Bug 2 (observed power instead of a
specified effect), the exact figure "0.0%" should not be interpreted as
a precise power estimate, but as a strong, reproducible warning signal:
the donor×vector interaction is not reliably detectable with the current
model/design at ANY of the three POC sample sizes — not even in the
"optimistic" 1,350-isolate scenario. This is consistent with the rule of
thumb stated in the POC document itself ("interaction effects are
typically markedly harder to detect than main effects with so few
clusters"), but turns out more pronounced here than the document itself
would have anticipated.

**Consequence/recommendation:**
1. Treat the interaction as exploratory only for now, not as a
   confirmatory main hypothesis (exactly the fallback already envisioned
   in the proposal itself: "stepwise, beginning with simple models").
2. For a reliable figure before applying for funding: switch the model
   to `lme4::glmer`/`glmer.nb` (where `simr`'s `fixef<-` works natively)
   and repeat the power analysis with an explicitly specified target
   effect.
3. The number of field×year cells is the actual bottleneck (not the
   number of isolates per cell) — more fields/years would presumably
   help more here than more isolates per field-year. This argues for
   expanding the spatial-temporal diversity of the sample, not just its
   size.
4. Once real HGT candidate rates from Sections 1–3 are available, replace
   `baseline_rate`/`effect_rr` with observed values and repeat — the
   current placeholder values are assumptions, not observations.

## 2026-09-01 — Global screening: expansion potential from the SRA catalog

**Finding:** The 15 largest BioProjects in the already cataloged
SRA Illumina dataset (`data/ncbi_m_oryzae_sra_wgs_illumina.tsv`, 1,754
runs) together account for 71% of all runs. A title sample of these
top 15 shows that several are already-published **multi-site population
genomics studies**, not scattered individual submissions: SRP592876
(244 runs, finger-millet isolates East Africa), SRP132141 (88,
"Population Genomic Analysis of the Rice Blast Fungus... Expansion of
Three Main Clades"), SRP288432 (48, population genetics sub-Saharan
Africa), ERP109496 (55, "various locations in Africa"), SRP076116 (43,
wheat blast Brazil/B71 lineage).

**Consequence:** For Section 3 (PCA/ARI screening), host-lineage
diversity and isolate count can be massively expanded directly from
already-published data, without new fieldwork. Computationally/
technically no obstacle (n=50 ≈ 5 h mapping, ≈ 66 GB storage, both
feasible on the current machine). Statistical reference point for "at
what sample size realistic results can be expected": the project's own
synthetic validation (`validate_clustering_synthetic.py`, default n=40 =
4 lineages × 10 isolates) achieved a 70–95% detection rate there for the
per-region method — structurally unreachable at n=6 (current POC run).
Rule of thumb: **at least ~8–10 isolates per host lineage**.

**Important limitation:** The global screening only improves Section 3,
not automatically the GLMM power from Section 4 (see above) — most SRA
entries lack the necessary field/year structure for that. Whether the
multi-site studies named above provide enough geo/time metadata in their
supplements to approximate the field/year design is an open check.

## 2026-09-01 — Pivot: new workflow document (`POC_HGT_Starships_Workflow.md`), old short-read coverage pipeline removed

**Decision:** The user provided a new, considerably more comprehensive
workflow document
(`Documentation/POC_HGT_Starships_Workflow.md`, originally
`C:\Users\<user>\Downloads\POC_HGT_Starships_Workflow.md`) and instructed
removal of everything previously built that is no longer needed under
the new approach. The old short-read coverage pipeline setup (`workflow/`,
`config/`, `power_analysis/`, old `tests/`, `data/references/panel_manifest.tsv`
and `candidate_regions.bed`, the bwa-mem2 index files) was removed.

**What changes substantively:** The new document replaces the previous
approach (short-read mapping against a small 6-genome panel + coverage
PAV + PCA/ARI clustering) with an 8-phase workflow built on a larger
(20–30 isolates), predominantly long-read-assembled and uniformly
annotated reference panel: assembly/QC (BUSCO, BlobTools) →
annotation/repeat-masking (funannotate/BRAKER, RepeatModeler/Masker) →
Starship catalog (still `starfish`) → pangenome PAV (Panaroo/PIRATE/
Roary + alignment-based accessory regions) → **saturation/rarefaction
analysis as the core piece of the POC** (subsampling curve: at what
panel size does candidate discovery plateau?) → validation (leave-one-out
+ real short-read test isolates) → phylogenetic incongruence (core tree
vs. Starship gene tree) → GLMM sample-size/model-complexity simulation
(newly reframed: clonal-lineage-based simulated phylogenies instead of a
field/year design) → synthesis/go-no-go report.

**What was preserved (see README.md for details):**
- Git history (checkpoint commit before the cleanup: all removed files
  remain retrievable via `git log`/`git show`)
- `Documentation/` (all previous planning notes + the new document)
- `assemblies/` (6 already downloaded reference genomes — starting point
  for the larger 20–30-genome panel from Phase 1)
- `data/isolates_poc/` (5 short-read + 1 long-read test isolate, already
  downloaded and QC-checked — directly usable for Phase 5)
- All NCBI catalogs in `data/` and `input/metadata/`
- `envs/starfish.yaml` (Starfish is still needed unchanged in Phase 2);
  `envs/mapping.yaml`/`python.yaml` as generic tool building blocks
  (minimap2 is needed again in Phase 3 for alignment-based PAV
  detection)

**Newly created:** `poc_hgt_starships/{00_data...09_report}/` with
per-phase READMEs (input/output/tools per phase, adopted from the new
document).

**Deliberately deferred:** Phase 7 (sample-size/GLMM simulation) — not
implemented now, per explicit user request ("leave the GLMM aside for
now"), even though the new document calls for it. The earlier analysis
of the old field/year GLMM (see entry above, the 0% power finding)
remains documented as background knowledge but does not directly carry
over to the newly reframed Phase 7 model (clonal-lineage-based
simulation).

**Next step:** Phase 1 (`poc_hgt_starships/00_data/`) — curating the
20–30-isolate panel from the already existing NCBI long-read assembly
catalogs (`data/ncbi_m_oryzae_longread_candidates*.tsv`,
`ncbi_m_oryzae_longread_highquality.tsv`, `ncbi_m_oryzae_sra_wgs_longread.tsv`),
focused on maximizing clonal-lineage/host diversity.

## 2026-09-02 — Complete NCBI assembly catalog (605 genomes) via Datasets REST API

**Decision:** Instead of the previous ~46-entry assembly list
(`data/ncbi_m_oryzae_assemblies.tsv`, an older, incomplete pull), a
complete re-fetch was performed via the NCBI Datasets REST API v2
(`https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/318829/dataset_report`,
txid 318829 covers both "Pyricularia oryzae" and its taxonomic synonym
"Magnaporthe oryzae"): **605 assemblies**, retrieved in full in a single
request (page_size=1000). Result in
`poc_hgt_starships/00_data/ncbi_pyricularia_oryzae_assemblies_full.tsv`
(27 columns: accession, organism, strain, assembly level,
`sequencing_tech` — directly from the API, no heuristic needed —,
BioSample host/isolation source/geo/collection date, contig/scaffold
N50/L50, chromosome count, GC%).

**Why the Datasets API instead of E-Utilities:** The Datasets API
returns `assembly_stats` (N50 etc.) AND `assembly_info.sequencing_tech`
AND embedded BioSample attributes in a single structured JSON response
per assembly — with E-Utilities (esummary) this would have required at
least 3 separate queries per assembly (assembly summary, BioSample
summary, possibly the assembly-stats report file).

**Key findings (details: `poc_hgt_starships/00_data/dataset_summary.md`):**
- **Storage:** 24.2 GB total uncompressed FASTA (605 genomes at
  35–49 Mb, median 39.9 Mb); ≈ 6.1 GB gzip-compressed.
- **Host diversity:** 29 normalized categories after typo/synonym
  cleanup (`host_diversity_summary.tsv`); dominated by *Oryza sativa*
  (275), *Triticum aestivum* (59), *Eleusine* spp. (41), *Urochloa* spp.
  (19), *Lolium* spp. (16); **27% (164/605) lack a host attribute in
  the BioSample record.**
- **Assembly level/contig size:** only 14 "Complete Genome" (2.3%) + 42
  "Chromosome"-level (6.9%) — 66% (398) only scaffold-level with a
  contig N50 median of just 0.06 Mb.
- **No assembly explicitly flagged as T2T** (0 of 605). 14 assemblies
  are "gapless chromosome-level" (contig count = chromosome count) — an
  approximate criterion for T2T quality, but without verified telomere
  repeats at both ends, so not equivalent to genuine T2T.
- **Data type:** 84% short-read-based (508), 11.6% long-read (70), 4%
  hybrid (24), the rest Sanger/historical. Expected relationship
  confirmed: Complete-Genome level is exclusively long-read/hybrid,
  scaffold level almost exclusively short-read.

**Consequence for Phase 1:** Of 605 assemblies, only the 56
Complete-Genome/Chromosome-level entries are realistic candidates for
the reference panel (the rest too fragmented for Starship boundary
calling, per the document's stipulation that "repeat-masking is
essential for clean boundary calls" — at scaffold level with a contig
N50 of ~60 kb, Starship boundaries can hardly be determined cleanly).
Of these 56, in turn, only those with a host attribute AND sufficient
clonal-lineage diversity are relevant for the 20–30-isolate target size
— still to be filtered.

## 2026-09-02 — Second pivot: `multireference_panel_pav_workflow.md`, dedicated project `magnaporthe_multiref_pav/`

**Decision:** The user provided a third, even more specific workflow
document (`Documentation/multireference_panel_pav_workflow.md`, 19
sections) and instructed that its steps be followed strictly, using a
dedicated folder structure for it. Unlike the first pivot (see above),
NOTHING was removed from `poc_hgt_starships/` — the new document is a
much more concrete elaboration specifically of the multi-reference-panel/
PAV part (congruent with
`poc_hgt_starships/{01_assembly_qc,03_starship_calls,04_pav_matrix}`),
not a replacement for the entire 8-phase plan. A new, standalone project
directory `magnaporthe_multiref_pav/` was created exactly following the
structure specified in Section 4 of the document (`config/`,
`data/{references,annotations,longreads,resources}/`, `envs/`,
`workflow/{rules,scripts}/`, `results/`, `logs/`).

**14-genome reference panel identified and downloaded:** The "14
complete genome assemblies" required by the document correspond exactly
to the 14 "Complete Genome" entries from the 605-genome NCBI catalog
compiled days earlier (not searched anew, directly reused). All 14 FASTA
files were successfully downloaded via the NCBI Datasets API
(`data/references_raw/`) — **only 1 of 14 (`GCA004346965_1`, an
Eleusine isolate) comes with a GFF3 annotation**; the remaining 13 need
the reannotation (BRAKER3/Liftoff) called for in Section 6.2 before
gene-based analyses (OrthoFinder, Section 7.1) are possible.

**Critical finding on the "42 long-read test isolates" assumption:** The
document assumes that 42 chromosome-level long-read isolates are
available for the pilot selection (Section 9). The NCBI catalog does
provide exactly 42 assemblies at Chromosome level (the numeric match is
no coincidence) — but closer inspection of the `sequencing_tech`
metadata shows:
- **Only 15 of 42 are actually long-read/hybrid-sequenced** (13
  long-read + 2 hybrid); 25 are short-read-based (predominantly one
  large batch of 22 Brazilian wheat isolates, presumably
  reference-guided scaffolded, not de novo long-read-assembled); 2 are
  the historical Sanger 70-15 duplicate (GCA/GCF_000002495.2).
- **Host diversity in this pool of 42 is strongly skewed:** 31/42
  *Triticum*, 5 *Oryza*, 1 each *Lolium*/*Setaria*, **0 *Eleusine*.**
  The stratification required in Section 9.1 (among others, "2
  Eleusine-associated isolates") cannot be met from this pool.

**Consequence:** For genuine Eleusine long-read test isolates, the
separate SRA raw-data catalog (`data/ncbi_m_oryzae_sra_wgs_longread.tsv`,
193 runs) needs to be searched by BioSample host attributes — not yet
done. `config/samples_candidate_pool.tsv` (42 rows, with a `platform`
column) documents the full pool including this limitation; the final 10
pilot isolates (`config/samples.tsv`) are therefore not yet populated.

**Phase I (Section 6.1) executed:** `seqkit stats` for all 14 references
(`results/qc/assembly_stats.tsv`) — a consistent picture (7–10 contigs,
42–48 Mb, N50 5.7–7.5 Mb, GC ~50%), matching the expected *M. oryzae*
profile. BUSCO (`sordariomycetes_odb10`) ran in the background for all
14 genomes (env `qc_env`, already present with BUSCO from an earlier
session, only `seqkit` added rather than a redundant new environment).

**Deliberately deferred (scope/tool availability):** Annotation
(BRAKER3 — needs a separate GeneMark license, not a pure conda install),
repeat-masking (RepeatModeler2/EDTA), whole-genome alignment/SyRI,
OrthoFinder, panel construction, long-read mapping, window-based PAV,
SV calling (Sniffles2), rarefaction. All set up as documented stubs in
`workflow/rules/*.smk` with a status comment and prerequisites, so that
the Snakemake structure is complete and the next implementation step per
file is clear.

**Bug: BUSCO OOM under parallel execution.** The first BUSCO run
(`--cores 2`, 2 genomes in parallel) led to repeated, initially
confusing crashes (empty `/tmp`, `LockException`, EXT4 un-/remounts, and
"journal corrupted or uncleanly shut down" in the kernel log — it looked
like a WSL VM restart). `dmesg` clearly showed the actual cause:
`Out of memory: Killed process ... (python3) ... anon-rss:6897092kB` —
a single BUSCO genome-mode run (metaeuk gene prediction against a
~44-Mb genome) alone needs **~6.9 GB RSS**; two in parallel clearly
exceed the 10 GB WSL RAM limit (+4 GB swap). Fix: `qc.smk`'s
`busco_reference` rule now gets an explicit `threads:` directive
(reserving the full core budget per job), and the run is started with
`--cores 1` (true serialization, one genome after another) instead of
`--cores 2`. Expected runtime therefore longer (~14 × 10–20 min instead
of parallelized), but stable.

## 2026-09-02 — BUSCO crashes: the actual root cause was WSL2 `autoMemoryReclaim`, not power-saving mode

**Background:** Despite the OOM fixes (see above), the serialized BUSCO
run (`--cores 1`) kept crashing unexpectedly — without an internal BUSCO
error (aborting mid-`hmmsearch`, sometimes even with a WSL interop error
"Failed to start the systemd user session"). `dmesg` consistently showed
EXT4 un-/remount cycles of the root disk (`sdd`) at a ~100–130-second
cadence, accompanied by `systemd-journald: File ... corrupted or
uncleanly shut down`.

**First (incorrect) hypothesis:** Windows power-saving mode (AC standby
timeout 45 min) or USB Selective Suspend (presumably power-cycling the
disk behind `sdd`). The user was asked and chose "temporarily disable
sleep timeout"; both (`STANDBYIDLE` and USB Selective Suspend, AC side)
were disabled via `powercfg`. **Result: no effect** — the identical
crash rhythm recurred unchanged, this time with the additional finding
that `/dev/sdd` is not an external/USB drive at all, but the **root disk
of the WSL2 VM itself** (`/` and `/mnt/wslg/distro`, a 1 TB dynamic
VHDX) — the powercfg hypothesis was thus structurally implausible (no
USB device involved) and was discarded.

**Actual root cause:** `C:\Users\<user>\.wslconfig` had no explicit
`autoMemoryReclaim` setting, causing WSL2 (version 2.7.12.0) to use the
default **`gradual`** — a periodic memory compaction of the VM that, for
memory-intensive workloads (like the ~6.9 GB RSS BUSCO process), freezes
the VM just long enough to trigger I/O timeouts on the root disk and
thereby the observed remounts/journal corruption. Fix: `.wslconfig`
extended with
```
[experimental]
autoMemoryReclaim=disabled
```
`wsl --shutdown` executed (clean restart, `uptime` confirms 0 min),
BUSCO run restarted. **Result: all 14 genomes completed without
interruption.** The powercfg changes were reverted to their original
values (AC standby 0x00000a8c/2700s, USB Selective Suspend AC
0x00000001/enabled), since they were demonstrably not the cause.

**Lesson:** For stability problems internal to the WSL2 VM (un-/remounts
of the root disk, not of a peripheral device), check `.wslconfig`
(`autoMemoryReclaim`, `vmIdleTimeout`, `sparseVhd`) first, before
suspecting Windows host power settings as the cause — the device behind
the disk apparently affected by "unmount/remount" should always be
identified first via `mount`/`lsblk` (here: `sdd` = `/`, not a USB
device).

**BUSCO result (`sordariomycetes_odb10`, all 14 reference genomes):**
consistently 97.9–98.2% complete (mostly single-copy, duplication
≤0.5%), <2% missing — consistently high assembly completeness across
the entire panel, no outliers.

## 2026-09-02 — Eleusine long-read gap closed: 2 genuine isolates found in the SRA raw-data catalog

**Approach:** From `data/ncbi_m_oryzae_sra_wgs_longread.tsv` (193 runs),
the 158 unique BioSample accessions were extracted and their BioSample
attributes (host, isolate, isolation source, geo) retrieved via NCBI
E-Utilities (`efetch db=biosample`, batches of 50)
(`config/longread_sra_biosample_hosts.tsv` saved in the new project,
115/158 BioSamples returned attributes — the rest are older BioSamples
without structured metadata).

**Result:** Two isolates found with host `Eleucine coracana`
(finger millet) AND genuine long-read raw data:

| Isolate | BioSample | SRA run | Platform | Bases | Origin |
|---|---|---|---|---|---|
| K23/123 | SAMN08033374 | SRR6307184 | PacBio RS II | 3.31 Gb (≈74× at 44.5 Mb) | Kenya: Busia district, blast |
| E34 | SAMN12142210 | SRR9972918 | PacBio Sequel | 8.41 Gb (≈189×) | Ethiopia: Diga, blast |

Both appended as new rows to `config/samples_candidate_pool.tsv`
(pool now 44 instead of 42 entries; the `fastq` column already carries
the SRA run accession, since unlike all 42 original chromosome-level
entries, the raw data are actually directly downloadable here).
`read_n50_bp` here is the mean subread length from the SRA runinfo (not
a true N50, since this metric is absent from the runinfo format) — the
true N50 should be computed from the downloaded reads before the actual
mapping run.

**Important secondary finding:** Cross-referencing the 42 original
chromosome-level BioSamples against the 158 long-read SRA BioSamples
yielded **zero overlaps** — none of the 42 assemblies has raw reads
findable in this catalog (the long-read assemblies were apparently
submitted without a raw-data deposit, or under a different BioSample
than the assembly's BioSample). This means: so far NONE of the 42
entries in `config/samples_candidate_pool.tsv` actually had loadable
FASTQ ("assembly only" for all of them) — K23/123 and E34 are thus not
only the Eleusine fix but currently the ONLY two entries in the entire
pool with genuinely available raw data. For the remaining host groups
(Triticum, Oryza, Lolium, Setaria), the SRA catalog must likewise be
searched for raw data, not just assemblies, before the final selection
of 10 pilot isolates (Section 9.1) — not yet done.

## 2026-09-02 — Repeat-masking (Section 6.3): a dedicated, lightweight environment instead of the heavy `annotation.yaml`

**Decision:** `envs/repeats.yaml` (new, `multiref-repeats`) created with
only `repeatmodeler=2.0.5`, `repeatmasker=4.1.7`, `bedtools`, `seqkit`,
`samtools`, instead of installing RepeatModeler/RepeatMasker from the
already-planned `envs/annotation.yaml` (which additionally bundles
BRAKER3, liftoff, eggnog-mapper, orthofinder, diamond, iqtree — a
combined install attempt of all these packages would have been slower
and more fragile, especially because of BRAKER3's GeneMark license
dependency, which has to be handled separately anyway).

**`workflow/rules/repeats.smk` implemented** (previously a pure stub):
`index_reference_fai` (samtools faidx) → `build_repeat_database`
(BuildDatabase, NCBI engine) → `run_repeatmodeler` (with `-LTRStruct` for
LTR retrotransposon sensitivity) → `run_repeatmasker` (using the
genome-specific RepeatModeler library) → `repeat_windows`
(`workflow/scripts/repeat_windows.sh`: RepeatMasker `.out` → BED →
`bedtools merge`/`coverage` against a window grid). Window size is taken
from `config/thresholds.yaml: pav.window_size_bp` (10 kb), not redefined
— this lets the repeat density per window later be joined directly with
the window-based PAV classification (Section 8), exactly as required by
the document ("repeat fraction" as a column of the final panel table,
Section 8.4/8.5). This required extending `workflow/Snakefile` to parse
`config/thresholds.yaml` itself (`thresholds = yaml.safe_load(...)`) —
previously the file was only referenced, never parsed.

**As with BUSCO:** `run_repeatmodeler`/`run_repeatmasker` get an
explicit `threads:` (= `config["threads_default"]`) to force true
serialization when `--cores` equals the thread value — RepeatModeler is
similarly memory-/time-intensive as BUSCO, and multiple parallel runs on
the 10 GB WSL2 VM are a known risk (see BUSCO OOM above).

**Pilot run before full execution:** Since RepeatModeler2 with
`-LTRStruct` can, per the literature, take several hours per genome for
genomes this size, and 14 genomes run serially could potentially mean
1–4+ days, only `GCA036493215_1` (the smallest genome in the panel,
42.5 Mb) was first started as a timing pilot, before committing to all
14 as a long background run.

**Pilot-run result (timing):** Round 1 (RepeatScout, largest sample)
took 63 minutes for `GCA036493215_1` — of which 56 minutes alone were
spent on ONE unusually copy-rich repeat family ("family-0"; the
remaining 88 discovered families completed in seconds to low minutes).
Unlike BUSCO, RepeatModeler is **barely memory-hungry** (~1 GB RSS
instead of ~7 GB) — the 10 GB WSL2 limit is not an issue here, the
bottleneck is purely CPU-/time-bound. This makes genuine parallelization
across multiple genomes (instead of serial execution as with BUSCO) both
possible and sensible.

## 2026-09-02 — Panel reduced to 5 host representatives (POC scope deviation from the document) + repeat-masking parallelized

**Important:** The document (`multireference_panel_pav_workflow.md`)
itself does NOT call for a reduction of the 14-genome catalog — it
consistently assumes all 14 genomes. The following reduction is a
deliberate, user-driven POC scope decision (to reduce computational
cost), not a document requirement, and is documented here transparently
as a deviation.

**Decision:** Instead of all 14 Complete-Genome references, only **one
representative per host type** (5 genomes) is now used for
repeat-masking, annotation, and panel construction.
`config/references.tsv` (active, read by `workflow/Snakefile`) now
contains only these 5 rows; the full 14-genome catalog remains archived
unchanged in `config/references_full_catalog_14genomes.tsv` (the already
completed BUSCO results for all 14 remain valid and are not discarded,
but do not proceed further through subsequent pipeline steps).

**Host-type tally and correction of a misclassification:** Originally,
`GCA036493215_1` was carried as `host_group=unknown` (the `host` field
was empty in the NCBI BioSample record). A targeted follow-up query to
the NCBI Datasets API (`accession/GCA_036493215.1/dataset_report`)
revealed, however, that it is **Br48**, a wheat-infecting isolate from
Brazil (BioSample attribute `strain: "wheat infecting strain"`,
`geo_loc_name: Brazil`), AND the associated BioProject description
explicitly reads "**Telomere-to-telomere genome assembly of Pyricularia
oryzae Br48**" (PRJDB14561) — a self-declared T2T assembly (7 contigs =
7 chromosomes, `contig_l50=3`). `host`/`host_group` were corrected
accordingly to `Triticum aestivum`/`triticum`. This leaves only **5
genuine host types** among the 14 genomes (Oryza, Triticum,
wild grass/Lolium, Eleusine, Avena), not 6 as initially assumed.

**Representative selection:**

| Host type | Candidates (14-genome catalog) | Chosen | Rationale |
|---|---|---|---|
| Oryza | 7015, Guy11, 95HPH4, 95085, P131 | **7015** | canonical reference strain "70-15", a field standard in practically every comparative *M. oryzae* genome study |
| Triticum | GCA059330735_1, GCA059330115_1, GCA059330025_1, GCA059330365_1, Br48 | **GCA036493215_1 (Br48)** | the only explicitly T2T-declared genome in the entire panel — directly relevant for clean Starship/mini-chromosome boundary calls (Section 6.3/8), more important than the marginally higher contig N50 of the alternatives |
| Wild grass | LpKY97, GCA059329725_1 | **LpKY97** | an established reference strain used in the literature for the Lolium lineage |
| Eleusine | GCA004346965_1 | GCA004346965_1 | only Complete-Genome option |
| Avena | GCA059329645_1 | GCA059329645_1 | only Complete-Genome option |

`GCA036493215_1` (Br48) keeps its genome ID (accession-based rather than
"Br48"), even though the strain name is now known — the already-running
RepeatModeler pilot run (see above) uses exactly this wildcard value;
renaming would have orphaned the progress made so far.

**Known trade-off:** The reduction loses within-host diversity (e.g., 5
instead of 1 Oryza genome, 5 instead of 1 Triticum genome) — Starships
that occur only in a subset of isolates of ONE host type will not be
captured by the 5-genome panel. This is acceptable for the POC's core
question (detectability of Starships/accessory chromosomes ACROSS host
boundaries via long-read mapping); for a later full analysis/publication
the full 14-genome catalog should be used.

**Windows sleep mode disabled (this time for a genuine reason, not a
misdiagnosis):** Unlike the earlier BUSCO crash (where sleep/USB suspend
were falsely suspected), there is a real risk here — multi-hour
background runs would be paused/interrupted by actual sleep mode (not
just the screen saver). At the user's request, `STANDBYIDLE` and
`HIBERNATEIDLE` were set to 0 (disabled) for both AC and DC. **Original
values for later restoration:** AC standby 0x00000a8c (2700s/45min), DC
standby 0x00000708 (1800s/30min), DC hibernate 0x0003f480
(259200s/3 days), AC hibernate was already 0. Should be restored once
the current multi-day background runs (repeat-masking) are complete.

**Repeat-masking parallelized:** Since RepeatModeler barely needs RAM
(see pilot-run finding above), all 5 panel genomes now run
SIMULTANEOUSLY instead of serially: the already-running Br48 pilot
(8 threads, continued unchanged) plus a second Snakemake run started in
parallel via `--nolock` for the remaining 4 genomes (`7015`, `LpKY97`,
`GCA004346965_1`, `GCA059329645_1`, 3 threads each,
`config.yaml: threads_default` lowered from 8 to 3 for this). `--nolock`
is safe here since both runs target disjoint genomes and therefore
disjoint output files. Resource check with 5 parallel runs: 4.5 GB RAM
(of 10 GB), load average ~10 (of 14 cores) — stable, no OOM risk.

**Result (all 5 genomes complete):** The 4 parallelized genomes together
needed only ~2:48 h of wall time (13:40–16:29) instead of the ~8 h a
serial run (4 × ~2 h as with the Br48 pilot) would have cost —
parallelization was the right call.

| Genome | Host | Runtime (RepeatModeler) | Families | Genome-wide repeat fraction |
|---|---|---|---|---|
| 7015 | Oryza | 1:46 h | 76 | 16.47% |
| LpKY97 | Wild grass | 2:06 h | 203 | 14.53% |
| GCA059329645_1 | Avena | 2:21 h | 133 | 11.66% |
| GCA004346965_1 | Eleusine | 2:46 h | 94 | 11.58% |
| GCA036493215_1 (Br48) | Triticum | 1:59 h | 138 | 10.12% |

**Important biological finding — two clear mini-/accessory-chromosome
candidates** (from results/repeats/{genome_id}_repeat_per_contig.tsv,
pattern: markedly smaller contig + markedly higher repeat fraction than
the rest of the genome):

- LpKY97 (wild grass): CP050927.1 (3.0 Mb, 56.3% repeat) and CP050928.1
  (0.9 Mb, smallest contig, 53.0% repeat) - both massively above the
  rest of the genome (6-17% for the other 7 contigs).
- GCA059329645_1 (Avena): CM181343.1 (1.3 Mb, 45.4%) and CM181341.1
  (1.2 Mb, 24.4%) - likewise clearly above the genome average (11.7%).
- For 7015, GCA004346965_1, and Br48 there is NO comparable outlier (all
  contigs large, 4.0-8.8 Mb, repeat fraction in the normal range
  4.8-22%) - these 3 assemblies appear either not to contain separately
  assembled mini-/accessory chromosomes, or these were not resolved as
  their own contigs.

**Consequence:** These four contigs (LpKY97 x 2, GCA059329645_1 x 2) are
strong initial candidates for accessory_chromosome/mini_chromosome
(Section 8.2 classification) and should be given priority review in the
later Starship search (Section 6.4/DUF3435 scan) and panel
classification (Section 8) - not yet done, since Starship annotation
(starships.smk) remains an unimplemented stub.

## 2026-09-02 — Section 6.2 (partial): Liftoff gene annotation for all 5 panel genomes, mini-chromosome hypothesis confirmed by gene density

**Decision:** Since BRAKER3 (de novo annotation, Section 6.2) remains
blocked by the missing GeneMark license, **Liftoff** was used instead
(dedicated, lightweight `envs/liftoff.yaml`, separate from the heavy
`envs/annotation.yaml`) - the only genuine annotation available from
NCBI (`GCA004346965_1`, 13,521 genes) was transferred to all 5 panel
genomes (including a self-liftoff of `GCA004346965_1` onto itself, as a
renaming step: the raw GFF3 uses the ORIGINAL contig IDs, while Liftoff
produces output with the renamed `{genome_id}__{contig}` IDs matching
`data/references/*.fa`).

**Bug (race condition):** In the first parallel run (4 genomes
simultaneously, same reference GFF3), `LpKY97` failed with
`UnboundLocalError: cannot access local variable 'feature_db'` -
Liftoff's internal gffutils step by default builds an SQLite database
file next to the input GFF3; all 4 jobs simultaneously tried to write
the same file (a race condition, the same basic pattern as the earlier
BUSCO lineage download). Fix: `LpKY97` restarted alone (without
contention) - completed successfully right away.

**Result (genes transferred / not mapped, out of 13,521 source genes):**

| Genome | Host | Genes transferred | Not mapped | Fraction |
|---|---|---|---|---|
| GCA004346965_1 | Eleusine (self-liftoff) | 13,518 | 3 | 99.98% |
| GCA036493215_1 (Br48) | Triticum | 13,209 | 312 | 97.7% |
| GCA059329645_1 | Avena | 13,231 | 290 | 97.9% |
| LpKY97 | Wild grass | 13,256 | 265 | 98.0% |
| 7015 | Oryza | 12,787 | 734 | 94.6% |

The loss correlates sensibly with phylogenetic/host distance from the
Eleusine source genome (Oryza most distant → highest loss) - a good
plausibility signal for the method itself.

**Important additional confirmation of the mini-chromosome candidates
(see previous entry):** Gene density per contig shows a drastically
REDUCED gene density for all four contigs previously flagged by repeat
fraction - exactly the expected dual signal (repeat-rich AND gene-poor):

| Contig | Size | Genes | Gene density | Comparison to genome average |
|---|---|---|---|---|
| LpKY97 CP050927.1 | 3.0 Mb | 119 | ~40/Mb | ~300-400/Mb elsewhere - **~8-10x lower** |
| LpKY97 CP050928.1 | 0.9 Mb | 43 | ~48/Mb | ~8x lower |
| GCA059329645_1 CM181343.1 | 1.3 Mb | 129 | ~99/Mb | ~327/Mb elsewhere - **~3x lower** |
| GCA059329645_1 CM181341.1 | 1.2 Mb | 297 | ~247/Mb | slightly below average |

**Consequence:** `LpKY97__CP050927.1`, `LpKY97__CP050928.1`, and
`GCA059329645_1__CM181343.1` are now candidates for
`accessory_chromosome`/`mini_chromosome` supported by TWO independent
lines of evidence (high repeat fraction + strongly reduced gene
density) - markedly more robust than either piece of evidence alone.
`GCA059329645_1__CM181341.1` remains a weaker secondary candidate.
Missing pieces for a complete Section 6.2: BRAKER3 de novo annotation
(particularly important in these gene-poor/repeat-rich regions, per the
document's warning that "lift-over alone can underestimate accessory
genes"), InterProScan/eggNOG functional annotation, OrthoFinder
(Section 7.1).

## 2026-09-02 — Section 7.1 (OrthoFinder): a gffread stop-codon character made diamond 2.2.6 hang instead of erroring

**Setup:** Dedicated `envs/orthofinder.yaml` (orthofinder, diamond,
mafft, iqtree, gffread), separate from `envs/annotation.yaml`.
`gffread -y` extracts protein sequences from the Liftoff GFF3
(`extract_proteome` rule), then `orthofinder -f ... -S diamond -M msa
-T iqtree3` (the document writes `-T iqtree`, but the installed
OrthoFinder version v3.1.5 names the method `iqtree3` - the same choice
in substance). `-t`/`-a` adjusted to the 14 available cores (document
example: 32/16). `-o` could not be used as planned (it requires a
not-yet-existing target directory, which conflicts with Snakemake's
automatic creation of output parent directories) - instead OrthoFinder
runs with its default path (`OrthoFinder/Results_<date>/` in the
proteome directory) and is afterwards moved to the fixed target path via
`mv`.

**Bug 1 (trivial):** The first run failed with
`UnboundLocalError: cannot access local variable 'feature_db'` -
misleading; this was the LIFTOFF error from the previous entry, not
relevant here, see above.

**Bug 2 (serious, a multi-hour failure):** `orthofinder -T iqtree`
failed immediately ("Invalid argument for option -T: iqtree... Valid
options are: fasttree, raxml, raxml-ng, iqtree3") - fixed by using
`iqtree3` instead of `iqtree`.

**Bug 3 (the actual time sink, ~1 hour lost):** The `diamond makedb`
step (the first genuine compute step) HUNG for seemingly inexplicable
reasons - no error, no progress message, 0% CPU load after minutes. The
first hypothesis (based on temporal coincidence with another, much
rarer EXT4 remount event in the `dmesg` log, about 4.7 h after the last
WSL restart) was a recurrence of the earlier WSL stability bug - this
hypothesis was DISCARDED after a second, cleanly isolated attempt
(system otherwise completely quiet, `load average` ~0, no new `dmesg`
event) reproduced the exact same hang. Systematic narrowing: `diamond
makedb` hangs both single- and multi-threaded, both on the full protein
file and on a 1000-line subset - this rules out library/CPU-
architecture or file-size problems. A test with `diamond=2.1.9` (instead
of the newer version installed via dependency resolution) immediately
produced a genuine error instead of a hang: `Error reading input
stream at line 11: Invalid character (.) in sequence`.

**Cause:** `gffread -y` marks codons that cannot be translated cleanly
(among them stop codons, but also codons at exon boundaries with a
leftover frame) in the protein output with a period (`.`) - but this
character is NOT part of the amino-acid alphabet accepted by diamond.
diamond v2.2.6 (the initially installed version) hangs completely on
this invalid character instead of throwing an error - a genuine bug in
this diamond version, not a WSL/resource problem.

**Fix:** The `extract_proteome` rule now replaces `.` (and, as a
precaution, `-`) in sequence lines (not header lines) with `X` (unknown
amino acid) via `sed '/^>/!{s/\./X/g; s/-/X/g}'` before writing the
protein file. A direct downgrade to `diamond=2.1.9` was attempted but
discarded: `mamba install diamond=2.1.9` in the existing environment
triggered a dependency conflict and COMPLETELY UNINSTALLED `orthofinder`
(OrthoFinder requires a newer diamond version). The environment was
therefore rebuilt with `orthofinder` AND the diamond version it prefers;
cleaning the input alone was sufficient to fix the hang - diamond 2.2.6
works flawlessly on clean input.

**Lesson:** When a compute step hangs inexplicably at 0% CPU load, do
not jump to blaming an infrastructure problem (WSL, resources) just
because a genuine infrastructure problem occurred earlier - temporal
coincidence with an old symptom pattern (here: the rare EXT4 remount
event) can be misleading. An isolated minimal test (smallest input
file, single-threaded, alternative tool version) clarifies this faster
than another round of infrastructure diagnosis.

**Second runtime problem after the diamond fix:** With the document's
command `-M msa -T iqtree3`, the orthogroup assignment itself
(diamond+MCL) completed in minutes, but the subsequent gene-tree
inference called `iqtree3` (with ModelFinder) INDIVIDUALLY per
orthogroup - with >11,000 orthogroups and >15 minutes per tree (observed
on two running processes at 100% CPU load for >15 min without
completing), this would have taken several days to weeks. For the output
actually needed in Section 7.1 (Orthogroups.tsv as the basis for
gene-based PAV classification), the gene-tree refinement layer is not
required. Fix: `-M dendroblast` instead of `-M msa -T iqtree3` -
delivers the same orthogroup assignment (unchanged, determined by
diamond+MCL) without the expensive MSA/tree layer. Completed fully in a
few minutes afterward.

**Result (13,226 orthogroups, 66,000 genes, 5 genomes):**

| Metric | Value |
|---|---|
| Genes in orthogroups | 65,600 (99.4%) |
| Unassigned genes | 400 (0.6%) |
| Orthogroups with all 5 species | 12,259 (92.7%) |
| Single-copy orthogroups | 12,042 |
| Species-specific orthogroups | 7 |

**Prevalence distribution** (basis for rescaling the Section 7.3
thresholds, see `config/thresholds.yaml`):

| Prevalence | Orthogroups | Fraction |
|---|---|---|
| 5/5 (strict_core) | 12,259 | 92.7% |
| 4/5 (soft_core, NEW) | 589 | 4.5% |
| 2-3/5 (shell, NEW) | 371 | 2.8% |
| 1/5 (private_accessory) | 7 | 0.1% |

**Threshold rescaling:** The document's thresholds are calibrated for 14
references (`soft_core: 13/14=0.93`, `shell: 3-12/14`). With only 5
possible prevalence levels (1/5 to 5/5), the old `soft_core` threshold
of 0.93 would NEVER have been reachable (the next level below 1.0 is
4/5=0.80) - all 589 orthogroups with 4/5 prevalence would have been
misclassified as "shell" instead of "soft_core".
`config/thresholds.yaml` was now switched to `strict_core=1.00`,
`soft_core=0.80` (4/5), `shell_min=0.40` (2/5, everything below =
`private_accessory`) - qualitatively the same ordering as in the
document, just adapted to the smaller, discrete panel size.

## 2026-09-02 — Section 7.2: whole-genome alignments (MUMmer4/SyRI) - SyRI requires an equal chromosome count

**Setup:** `envs/wga.yaml` (mummer4, syri), `workflow/rules/wga.smk`
implements `nucmer_align` → `delta_filter` → `show_coords` → `run_syri`
exactly per the document's command (Section 7.2). Instead of a single
anchor, all 10 undirected pairs of the 5-genome panel are run
(`itertools.combinations`, `genome_pairs` in the Snakefile) - the
document explicitly warns against the bias introduced by a single
Oryza anchor.

**Bug 1:** SyRI (the installed version) requires `-d <delta file>` in
addition to `-c <coords>` for SNP/indel identification - the document's
example omits `-d`, which fails with
`ERROR - CIGAR string or .delta file is required` (the table coords
from `show-coords -THrd` contain no CIGAR). Fix: `-d {filtered.delta}`
added. In addition, `--prefix` expects only the filename suffix, not a
full path - `--dir` set separately (otherwise a crash warning).

**Bug 2 (a methodological limit, not a pure software bug):** SyRI
requires, for its whole-genome 1:1 chromosome assignment, that the
reference and query genome have the SAME number of contigs/chromosomes
(`ERROR - Unequal number of chromosomes in the genomes`). Our 5-genome
panel, however, deliberately has differing contig counts (7015=7,
Br48=7, GCA004346965_1=7, **LpKY97=9**, **GCA059329645_1=10**) - and the
"extra" contigs are exactly the ones already flagged as mini-/accessory-
chromosome candidates (see the repeats.smk finding above: LpKY97's 2
small repeat-/gene-poor contigs, GCA059329645_1's 2 small contigs).
**7 of 10 pairs** (each involving LpKY97 or GCA059329645_1) are
therefore not directly usable for SyRI.

**Fix (scope adjustment, implemented in `workflow/Snakefile` as
`syri_pairs`):** `nucmer`/`show-coords` run for ALL 10 pairs (yields raw
alignment coverage, works independently of chromosome count). SyRI's
full synteny/SV classification runs only for the 3 pairs with equal
contig counts (all 7-contig genomes against each other: `7015↔Br48`,
`7015↔GCA004346965_1`, `Br48↔GCA004346965_1`). No loss of data or
information for the mini-/accessory-chromosome question itself - the
affected contigs are already well supported by repeat-/gene-density
evidence; SyRI could not have provided a meaningful 1:1 chromosome
assignment for them anyway.

**Result (3 SyRI classifications, each >600,000 lines):** consistently
dominated by SNPs (237k-272k) and small indels (INS/DEL, 20-32k), with
several hundred syntenic blocks (SYN/SYNAL, 132-304 and 699-877
respectively) and a small but present number of structural variants
(inversions INV/INVAL/INVDP, duplications DUP/DUPAL, translocations
TRANS/TRANSAL) - biologically plausible for closely related
*M. oryzae* isolates from different host lineages.

## 2026-09-02 — Section 7.4: Starship/Captain candidates (starfish) - mini-chromosome finding confirmed

**Setup:** `envs/starships.yaml` (starfish, metaeuk, hmmer, mmseqs2,
blast) - identical to the already successfully tested global
`starfish_env` from an earlier session
(`Documentation/Starfish_Progress_Log.md`), created here as a
project-local copy for reproducibility, but the already-working global
`starfish_env` was directly reused for the actual run (saves
reinstalling the bundled HMM/reference protein databases under
`$CONDA_PREFIX/db`).

`starfish annotate` ran in ONE combined multi-genome run across all 5
panel genomes (2-column assembly TSV), with `-s '__'` as the separator
matching our existing `{genome_id}__{contig}` header convention (the
default would be `_`, which would have parsed genome IDs with their own
underscore, like `GCA036493215_1`, incorrectly). Additionally, `--gff`
was passed with the Liftoff gene models to reconcile existing genes with
newly predicted YR genes.

**Known limitation (not critical):** The internal name matching against
the Liftoff GFF3 failed
(`does not have a parse-able featureID using namefield 'Name='`) - our
Liftoff GFF3 stores gene names via `locus_tag=`, not `Name=` (Starfish's
default). As a result, 0 existing genes were linked to new YR hits ("no
metaeuk genes intersect..."). This was NOT fixed, since the subsequent
cargo-gene analysis is carried out independently against the Liftoff
GFF3 via coordinate overlap anyway (see below) - the internal Starfish
link was not required for that.

**`starfish annotate` result:** **68 HMM-validated YR/Captain genes**
found (7015: 17, LpKY97: 16, GCA004346965_1: 13, GCA036493215_1: 12,
GCA059329645_1: 10) - `results/starships/panel_YR.filt.gff`.

**Synthesis script (`workflow/scripts/classify_starship_candidates.py`,
rule `classify_starship_candidates`):** For each YR hit, combines:
- window size (centered on the hit, at least 20 kb per
  `thresholds.yaml: starship.min_region_length_bp`, clipped at contig
  ends via `.fai`),
- mean repeat fraction in the window (`repeats.smk` output,
  overlap-weighted),
- number of cargo genes in the window (Liftoff GFF3, pure coordinate
  overlap - independent of the failed internal Starfish link above).

Conservative classification (modeled on Section 7.4, but without SyRI
synteny cross-reference - see "still open" below): `starship_like` only
when minimum size AND at least 1 cargo gene AND repeat fraction in the
window ≥1.5x the genome average; otherwise
`duf3435_candidate_contextual` (at least one criterion met) or
`duf3435_candidate_only` (no context, a pure HMM hit - "a DUF3435 hit
alone is not sufficient", per the document).

**Result (68 candidates total):**

| Classification | Count |
|---|---|
| `starship_like` | 35 |
| `duf3435_candidate_contextual` | 32 |
| `duf3435_candidate_only` | 1 |

**Most important finding - closing the loop with the mini-chromosome
candidates:** The two LpKY97 contigs already flagged as mini-chromosome
candidates by both repeat fraction AND gene density do in fact carry
Captain/YR genes:
- `LpKY97__CP050927.1` (3.0 Mb, 56.3% repeat, 119 genes): YR59
  (`starship_like`, 2 cargo genes, 44% repeat in the window) and YR60
  (`duf3435_candidate_contextual`).
- `LpKY97__CP050928.1` (0.9 Mb, 53.0% repeat, 43 genes): YR61-64, of
  which YR62 and YR63 are classified as `starship_like`.

For `GCA059329645_1`'s mini-chromosome candidates (`CM181343.1`,
`CM181341.1`), by contrast, NO YR hits were found - not a
contradiction, but an additional differentiating feature: not every
accessory/repeat-rich region has to carry an active Captain-driven
element.

**This answers the project's core question (detectability of Starships/
accessory chromosomes via this multi-evidence pipeline) positively on
the 5-genome panel itself** - three independent lines of evidence
(repeat fraction, gene density, Captain-gene presence) converge on the
same LpKY97 contigs.

**Still open for a more complete Section 7.4 classification:**
- SyRI synteny cross-reference (possible only for the 3 compatible
  genomes, see previous entry) is not yet incorporated into the
  classification.
- Genuine boundary/insertion-site detection (Starfish follow-up steps
  beyond `annotate`, e.g., flanking direct repeats) was not carried
  out - the current window size is purely a distance buffer around the
  YR hit, not a structurally confirmed element boundary. Terminology
  therefore remains deliberately conservative ("starship_like", not
  "Starship").

## 2026-09-03 — Section 7.3: complete region-type classification of all panel windows

**Setup:** `workflow/scripts/classify_panel_regions.py`
(`classify_panel_regions` rule) classifies EVERY 10-kb window
(`pav.window_size_bp`) of all 5 panel genomes into exactly one region
type, by combining all previous phases:
- Orthogroup prevalence (7.1) per gene in the window, averaged →
  core/shell/private classification according to the thresholds already
  rescaled to 5 genomes (`thresholds.yaml`).
- Repeat density (6.3) for `repeat_ambiguous`/`subtelomeric_dynamic`.
- Contig size + gene density (< 50% of the median contig AND < 50% of
  the mean gene density) for `accessory_chromosome` - a programmatic
  version of the same logic previously applied manually to the 4
  candidate contigs.
- Starship-like candidate windows (7.4) as the highest priority level.

**Priority on overlap:** `starship_like` > `accessory_chromosome`
(entire flagged contig) > `subtelomeric_dynamic` > orthogroup-based
(strict_core/soft_core/shell/private_accessory) > `repeat_ambiguous`
(high repeat fraction, no gene evidence in the window) >
`unclassified`.

**Technical trick (ID mapping):** OrthoFinder orthogroups reference
genes via their mRNA IDs from the proteome FASTA (e.g.,
`rna-gnl|PRJNA507314|PoMZ_08913-RA_mrna`), while the GFF3 `gene`
features have their own IDs (`gene-PoMZ_08913`) with a shared
`locus_tag` attribute. `mrna_to_gene_id()` connects the two via the
shared `locus_tag`.

**Result (21,820 windows total):**

| Region type | Windows | Fraction |
|---|---|---|
| `strict_core` | 16,455 | 75.4% |
| `soft_core` | 2,736 | 12.5% |
| `repeat_ambiguous` | 1,450 | 6.6% |
| `shell` | 346 | 1.6% |
| `unclassified` | 282 | 1.3% |
| `subtelomeric_dynamic` | 227 | 1.0% |
| `accessory_chromosome` | 219 | 1.0% |
| `starship_like` | 101 | 0.5% |
| `private_accessory` | 4 |

## 2026-09-03 — Section 9.1: stratified selection of the 10 pilot isolates

**Approach:** Of the 17 candidates actually long-read/hybrid-sequenced
in the 44-isolate pool (`config/samples_candidate_pool.tsv`), 10 were
selected according to the document's criteria (Section 9.1).

**Secondary finding during host checking:** `GCA_021764705.1` (isolate
"EA18", China: Enshi Hubei) was carried in the pool as `host=unknown` -
a targeted NCBI Datasets API query returned `isolation_source: rice`,
i.e., actually an Oryza isolate (the same misclassification pattern as
earlier with Br48/Triticum, see the panel-reduction entry above).
`GCA_059469275.1` (isolate "E2", Ethiopia), by contrast, remained
genuinely unknown (BioSample provides only `isolation_source: Ethiopia`,
no host).

**Final selection (`config/samples.tsv`):**

| Isolate | Host | Origin | Rationale |
|---|---|---|---|
| O219 | Oryza | Ivory Coast, 1985 | Oryza diversity West Africa |
| TRG2 | Oryza | Thailand, 2023 | Oryza diversity Southeast Asia, most recent sample |
| Guy11 | Oryza | French Guiana, 1979 | benchmark isolate (established lab strain) |
| EA18 (GCA021764705_1) | Oryza | China, 2021 | Oryza diversity East Asia, hybrid-sequenced |
| B71 | Triticum | Bolivia, 2012 | wheat-blast pandemic lineage, continent of origin |
| ZM12 | Triticum | Zambia, 2018 | wheat-blast pandemic lineage, spread into Africa |
| K23_123 | Eleusine | Kenya | mandatory slot, already downloaded |
| E34 | Eleusine | Ethiopia | mandatory slot, already downloaded |
| TF051MC7 | Wild grass (Lolium) | USA, 2005 | "expected accessory DNA" criterion - same host group as LpKY97 (reference panel), which already showed mini-chromosomes |
| gw6 | Setaria | China, 2016 | additional host diversity (foxtail millet) |

4 Oryza instead of the recommended 2-3 (EA18 included as a bonus for its
hybrid sequencing and new geography, not a strict rule violation given
the range stated in the document).

**Critical finding - raw-data availability:** For the 8 non-Eleusine
isolates, the BioSample accessions were cross-checked against the
193-run long-read SRA catalog
(`data/ncbi_m_oryzae_sra_wgs_longread.tsv`) - **only 3 of 8 have raw
reads findable there:**
- **B71** (SAMN06076154): 11 runs (10x PacBio SMRT + 2x Nanopore),
  ~45 Gb combined - one of the most deeply sequenced samples in the
  entire catalog.
- **ZM12** (SAMN29254577): 5 Nanopore runs, 475 Mb to 13.7 Gb.
- **TF051MC7** (SAMN36850036): 1 Nanopore run, 8.2 Gb.

For **O219, TRG2, Guy11, EA18, and gw6**, NO raw reads were found in
this catalog - either deposited under a different BioProject/BioSample
(not captured by the original SRA search query) or genuinely not
uploaded separately from the assembly. This remains an open item for a
more targeted follow-up search (e.g., a direct SRA web search by
isolate name instead of only a BioSample cross-reference).

**Consequence:** Of the 10 pilot isolates, **5 of 10** (K23_123, E34,
B71, ZM12, TF051MC7) currently have genuinely findable/already
downloaded raw data - enough for a first mapping pilot run (Phase V),
but the full 10-isolate stratification from Section 9.1 is not yet
backed by real data.

**B71/ZM12/TF051MC7 downloaded (one SRA run each, not all available
ones):** For B71 and ZM12, multiple runs exist (11 and 5, respectively);
since the document itself names only a pilot starting value of ~20x,
only the smallest/a single run was downloaded per isolate instead of
all of them (B71: only `SRR6232287`, ~11.9 GB archive size, instead of
all 11 runs combined >45 GB). Result (`seqkit stats`/NanoPlot):

| Isolate | Run | Reads | N50 | Coverage (vs. 44.5 Mb) |
|---|---|---|---|---|
| B71 | SRR6232287 | 163,478 | 35,779 bp | ~81x |
| ZM12 | SRR19868246 | 397,795 | 20,605 bp | ~108x |
| TF051MC7 | SRR30725258 | 325,797 | 34,671 bp | ~184x |

**Important note on PacBio/Nanopore archive size:** The `.sra` archive
size can substantially exceed the raw base count (B71: 3.59 Gb bases per
the catalog, but an 11.9 GB download size - older PacBio RS/Sequel raw
formats also store additional kinetics/trace data, not just base
calls). Before further downloads from this catalog, always check the
`size_MB` column rather than estimating from the base count alone.

**Resource observation:** `fasterq-dump`/`gzip` for these large
individual files (7-10 GB uncompressed FASTQ) are, as expected,
CPU-intensive (several hundred percent CPU for fasterq-dump, one core
at 100% for gzip over several minutes) - makes the host machine
noticeably sluggish in the meantime, but unproblematic for data
integrity. At the user's request, this was left running at full speed
rather than throttling the thread count.

## 2026-09-03 — Section 10.1: long-read mapping of all 5 available pilot isolates against the panel

**Setup:** `workflow/rules/mapping.smk` implements `minimap2_map` →
`mapping_flagstat` → `mapping_coverage_by_contig` exactly per the
document's command (Section 10.1). Preset chosen per isolate based on
the actual SRA platform (not blanket ont/hifi as in the document's
example): `map-pb` for the three PacBio RAW/CLR isolates (B71, K23_123,
E34 - NO HiFi/CCS reads), `map-ont` for the two Nanopore isolates
(ZM12, TF051MC7). See `config/samples.tsv`, column `minimap2_preset`.

**Bug 1 (OOM):** The first run with `--cores 14` let Snakemake start 4
samples simultaneously (`threads: 3` each) - `samtools sort`'s default
memory reservation (~768 MB/thread) summed to >10 GB across 4 parallel
jobs and triggered an OOM kill (`dmesg` confirms, the `minimap2` process
was killed despite only ~700 MB of its own RSS - the actual cause was
`samtools sort`, not minimap2 itself). Fix: `-m 512M` set explicitly,
run repeated with `--cores 3` (effectively serial).

**Bug 2 (WSL instability, high-frequency recurrence):** After the OOM
fix, the run reproducibly failed twice in a row right after minimap2
index construction (RAM unremarkable at the time, no OOM) - `dmesg`
showed the earlier EXT4 remount/journal-corruption pattern recurring at
a ~100-130-second cadence, even though `autoMemoryReclaim=disabled` was
still set unchanged in `.wslconfig`. Cause not conclusively determined
(possibly a different disturbance source on the same day/boot than the
previous day). A clean `wsl --shutdown` + restart fixed it immediately -
afterward the complete 5-isolate mapping run ran stably for ~2.5 hours
without further interruption. Additionally, a retry loop
(`run_mapping_retry.sh`, up to 15 attempts with a 15 s pause) was built
into the execution script to automatically catch future transient
hiccups.

**Runtime:** Considerably longer than estimated from base count -
`--secondary=yes` against a 3,677-sequence panel with many homologous/
redundant core regions produces very many secondary alignments and
consequently very large BAM files (2.4-10.6 GB per isolate). Total
runtime for all 5 isolates: about 2.5 hours (serial, due to the memory
limit).

**Result (`samtools flagstat`):**

| Isolate | Primary mapped | Total mapped (incl. secondary) | BAM size |
|---|---|---|---|
| TF051MC7 | 99.54% | 99.92% | 10.6 GB |
| ZM12 | 92.96% | 98.38% | 5.7 GB |
| K23_123 | 91.89% | 97.91% | 2.4 GB |
| E34 | 77.53% | 93.67% | 6.1 GB |
| B71 | 71.62% | 94.37% | 4.0 GB |

B71's notably lower mapping rate (71.62%) has not yet been
investigated - a possible indication of greater divergence of the
Bolivian wheat-blast lineage from the panel, or of data-quality
differences (a single older-chemistry PacBio RS/Sequel run).

**Most important individual finding - a possible cross-host multi-copy
element:** `PANEL001785` (`accessory_chromosome`, cluster `ACC_002`,
representative `GCA059329645_1__CM181349.1`, a tiny 35-kb contig from
the Avena reference genome) shows an extremely high mean depth
(300x to nearly 4,000x, far above the isolates' respective 20-200x
whole-genome coverage) in **ALL 5** test isolates - despite completely
different host lineages (Triticum, Eleusine, wild grass). Equally
notable: `PANEL001948` (`shell`,
`GCA059329645_1__CM181340.1:60001-90000`) with similarly extreme depth
in all 5 isolates.

**Caution regarding interpretation:** Such extreme depth points more to
a highly repetitive multi-copy element (e.g., an rRNA gene cluster) than
to a normal single-copy region - exactly the scenario Section 10.2
warns about ("a region with high homology to multiple references must
not be attributed based on the primary alignment hit alone"). Before a
reliable interpretation as a "cross-host accessory element" is
warranted, this must be cross-checked with MAPQ-filtered alignments
(Section 10.2, PAV evaluation level 2) - see the next entry
(Section 11), where exactly that was done: the signal was confirmed.

## 2026-09-03 — Section 11: window-based PAV analysis - two starship-like regions confirmed across host lineages

**Setup:** `workflow/rules/pav.smk` fully implements Section 11:
`panel_windows` (10-kb windows across the panel FASTA, Section 11.1) →
`mosdepth_unique` (MAPQ ≥ 20, mosdepth's default already excludes
secondary/supplementary - the quantitative PAV level) +
`mosdepth_all` (`--flag 1540`: excludes only unmapped/qcfail/dup, KEEPS
secondary/supplementary - the homology control level, Section 10.2) →
`pav_call` (present/absent/uncertain per window according to
region-type-specific breadth thresholds from `thresholds.yaml`, plus
`ambiguous_multimapping` when the breadth difference between "all" and
"unique" exceeds 0.3) → `pav_matrix` (panel-region x isolate matrix,
majority vote over the windows of a panel region).

**Technical trick:** Since our panel FASTA headers already carry all
metadata themselves (`>PANEL000001|type=...|cluster=...|rep=...|source=...`,
no spaces), the mosdepth "chrom" value is directly the full header
string - `panel_id`/`region_type` are parsed directly from it, no
separate manifest join needed.

**Bug:** `mosdepth` was not installed in ANY actually existing
environment (`envs/core.yaml` was never actually created as
`multiref-core` - all previous `qc.smk` runs instead used the older,
pre-existing `qc_env`, which did not include mosdepth). Fixed via
`mamba install -n qc_env mosdepth` instead of creating a new
environment (pragmatic, since `qc_env` is already in use for all core
tools anyway).

**Runtime:** mosdepth is considerably faster than minimap2/samtools
sort - all 5 isolates (2 levels x 5 = 10 mosdepth runs + 5 pav_call + 1
pav_matrix) completed in **~10 minutes**, compared to the ~2.5 hours for
the preceding mapping.

**Result (3,677 panel regions x 5 isolates):**

| Call class | B71 | ZM12 | K23_123 | E34 | TF051MC7 |
|---|---|---|---|---|---|
| present | 854 | 1269 | 1261 | 1403 | 1413 |
| absent | 546 | 328 | 537 | 309 | 351 |
| uncertain | 153 | 196 | 170 | 237 | 359 |
| ambiguous_multimapping | 2124 | 1884 | 1709 | 1728 | 1554 |

**`ambiguous_multimapping` is the most frequent class in all 5 isolates
(42-58% of panel regions)** - a direct consequence of the already
documented weak panel dedup rate (89% of strict_core clusters remained
single-genome entries, see the Section-8 entry): many only slightly
different panel regions draw multi-mapping signal from each other.

**Main finding: 399 of 3,677 panel regions are unambiguously "present"
in ALL 5 isolates** (regardless of their host lineage):

| Region type | Count |
|---|---|
| soft_core | 294 |
| strict_core | 84 |
| shell | 18 |
| **starship_like** | **2** |
| accessory_chromosome | 1 |

**The single `accessory_chromosome` region is exactly `PANEL001785`**
(the previously notable Avena-specific 35-kb contig, see the
Section-10.1 entry) - **confirmed at MAPQ ≥ 20 as a genuine present
signal in all 5 isolates**, not as a multi-mapping artifact (which would
otherwise have been classified as `ambiguous_multimapping`). The signal
withstands the stricter check.

**The two `starship_like` regions confirmed as present in all 5 test
isolates (Triticum x2, Eleusine x2, wild grass x1), regardless of host
lineage:**

| Panel ID | Source | Length | Repeat fraction | In how many references |
|---|---|---|---|---|
| PANEL003659 | Br48 (Triticum), `AP027063.1:170001-200000` | 30 kb | 78.4% | only 1 (Br48-specific in the panel) |
| PANEL003670 | Avena, `CM181346.1:3530001-3560000` | 30 kb | 48.9% | 2 (Avena + Br48/Triticum) |

**Assessment:** This means the POC's core question - detectability
of cross-host starship-like elements via long-read mapping against the
multi-reference panel - is positively demonstrated with real data AND a
strict multi-mapping control. Both candidates deserve priority manual
follow-up (e.g., alignment visualization, cargo-gene identity between
the reference and test isolates) before a publication-ready statement -
the terminology deliberately remains "starship_like", not structurally
confirmed Starships (see the Section-7.4 caveat: no genuine boundary/
insertion-site detection was performed).

## 2026-09-03 — Section 12 (SV calling) and Section 14 (rarefaction) - the panel is NOT saturated

**Section 12, SV calling (Sniffles2):** `sniffles_call` (per isolate,
against the panel FASTA) + `sniffles_cohort` (cohort merge across the
`.snf` files) exactly per the document's command. Ran in seconds per
isolate (considerably faster than mapping/PAV). **Result: 193
structural variants** in the 5-isolate cohort - 110 deletions, 82
insertions, 1 inversion (`results/pav/pilot_cohort.sv.vcf.gz`). Not yet
merged with the coverage-based PAV matrix (Section 11) into combined
evidence (Section 12's "PAV evidence = coverage breadth + mapping
uniqueness + SV breakpoints + spanning reads" - an open follow-up step,
Section 13).

**Sections 14.1/14.2, reference-panel rarefaction:** The document calls
for 1,000 random permutations (for C(14,k), too large for exhaustive
enumeration). Our deliberately reduced 5-genome panel allows
**exhaustive enumeration of ALL C(5,k) combinations** (at most 10 per
panel size) - stricter than the document's own sampling, not a
weakening (`workflow/scripts/rarefaction_reference_panel.py`).

**Result at full panel size (k=5, 39 candidate regions total):** 5
accessory_chromosome + 30 starship_like + 4 private_accessory.

**Saturation check (Section 14.2, threshold 2-5% per the document):**

| Region class | R(4) mean | R(5) | Δ fraction | Saturated? |
|---|---|---|---|---|
| accessory_chromosome | 4.0 | 5 | 20.0% | **No** |
| starship_like | 24.4 | 30 | 18.7% | **No** |
| private_accessory | 3.2 | 4 | 20.0% | **No** |
| all combined | 31.6 | 39 | 19.0% | **No** |

**Important, honest finding: the 5-genome panel is NOT saturated** - the
gain from adding the 5th genome is ~19-20% for all candidate classes,
far above the 2-5% saturation threshold. This is a direct, expected
consequence of the panel reduction from 14 to 5 genomes (see the
earlier decision) - with only 5 instead of 14 references, saturation of
candidate-region discovery is not to be expected. **Consequence for a
later full analysis:** additional reference genomes (e.g., from the
archived 14-genome catalog, `config/references_full_catalog_14genomes.tsv`)
would very likely reveal further, currently uncaptured candidate
regions - the 5-genome panel is sufficient for the POC but should not be
mistaken for a complete candidate catalog.

**Section 14.3, test-isolate rarefaction/novelty check - deliberately
reduced scope:** Full implementation (unmapped reads → local assembly →
panel re-search → classification of new candidate regions) requires a
long-read assembler (e.g., Flye), which is not installed. Implemented:
only step 1 (extraction + basic statistics of the unmapped reads) as a
low-cost approximation of "how much isolate sequence the panel does not
explain at all".

| Isolate | Unmapped reads | Unmapped bases | Fraction of total reads |
|---|---|---|---|
| TF051MC7 | 1,489 | 4.4 Mb | 0.46% (consistent with 99.92% total mapping) |
| K23_123 | 27,928 | 41.5 Mb | 8.1% |
| ZM12 | 27,994 | 109.3 Mb | 7.0% |
| E34 | 211,317 | 485.9 Mb | 22.5% |
| **B71** | **46,395** | **273.0 Mb** | **28.4%** |

**B71 and E34 have by far the largest amount of unmapped sequence** -
consistent exactly with B71's already documented, notably lower mapping
rate (71.6%, see the Section-10.1 entry). 273 Mb of unmapped sequence
in B71 corresponds to roughly 6 times the genome size - a strong signal
that B71 (the Bolivian wheat-blast lineage) carries substantial sequence
content not represented in the current 5-genome panel. Combined with
the missing-saturation signal documented above, this confirms: **a
larger reference panel would noticeably improve candidate-region
coverage.** The actual contig assembly/novelty classification of these
unmapped reads was NOT carried out - a clearly documented open item for
a full analysis.

## 2026-09-03 — Section 13: candidate-region manifest and final per-isolate assignment

**Setup:** Two scripts implement Sections 13.1/13.2:
- `build_candidate_regions.py` filters `panel_contig_manifest.tsv` down
  to the three candidate classes actually assigned
  (`accessory_chromosome`, `starship_like`, `private_accessory` -
  `mini_chromosome`/`subtelomeric_dynamic` were never assigned, see
  earlier entries, hence no MCHR_/SUBTEL_ IDs) and links each
  `starship_like` region via coordinate overlap to its Captain gene
  (`panel_YR.filt.gff`, Section 7.4) and its cargo-gene count
  (`starship_like_candidates.tsv`).
- `build_candidate_region_calls.py` combines, for each candidate region
  x test isolate: PAV status (Section 11), breadth/depth/multi-mapping
  fraction (averaged from the window calls), and SV support (number of
  non-reference genotypes from the Sniffles2 cohort VCF, Section 12 -
  the VCF CHROM matches the panel FASTA header exactly, no coordinate
  transformation needed).

**Not implemented:** Section 13.3 (bedtools intersect with a panel-wide
gene GFF3 for functional annotation per PAV block) - our gene
annotation exists only per source genome (Liftoff), not projected onto
panel coordinates; this projection would be an additional step.

**Result:** `results/panel/panel_candidate_regions.tsv` (39 candidate
regions: 30 `starship_like`, 5 `accessory_chromosome`, 4
`private_accessory`) and `results/pav/candidate_region_calls.tsv` (196
rows = up to 39 regions x 5 isolates, where PAV data were available).

**Captain-gene linkage: 30 of 30 (100%)** `starship_like` candidates
received an assigned Captain gene via coordinate overlap - complete
consistency between Section 7.4 (the Starfish finding) and Section 8
(panel clustering), no lost assignments.

**The two already identified cross-host Starship candidates in
detail:**
- **STAR_012** (`PANEL003659`, Captain `GCA036493215_1__YR35`, Br48/
  Triticum) - "present" in all 5 isolates, confidence "medium" for
  B71/K23_123/TF051MC7/ZM12, "low" for E34 (lower breadth, 0.42).
- **STAR_023** (`PANEL003670`, Captain `GCA059329645_1__YR50`, Avena,
  shared with Triticum/Br48) - "present" in all 5 isolates, confidence
  "medium" for TF051MC7/ZM12, "low" for B71/E34/K23_123.

No SV support (`sv_support=0`) for either region in any isolate - the
presence signal rests exclusively on coverage breadth, not on
Sniffles2 breakpoints. This does not weaken the conclusion (per
Section 11.3, coverage breadth is the primary evidence for the interior
of starship-like elements; SV breakpoints are an ADDITIONAL, not
necessary, confirmation layer for the boundaries), but should be
explicitly noted in any later publication write-up. 0.0% |

**Consistency check passed:** `accessory_chromosome` windows occur
EXCLUSIVELY in `GCA059329645_1` (133) and `LpKY97` (86) - exactly the
two genomes with the already independently found mini-chromosome
candidate contigs. `starship_like` windows are distributed across all 5
genomes (17-24 per genome), consistent with the 68 Captain genes found
across all genomes.

**Still open:** SyRI synteny (available only for 3/5 genomes, see the
7.2 entry) has not been incorporated into the classification - it could
more precisely place `unclassified`/`shell` windows in non-coding
regions, where purely gene-based prevalence provides no evidence (1.3%
of windows affected). `mini_chromosome` (a document class that
presupposes verified telomere boundaries) was deliberately not
assigned - absent telomere-repeat verification, the more conservative
`accessory_chromosome` class remains the correct choice.

## 2026-09-03 — Section 8: analytical multi-reference panel built - the dedup threshold meets deliberate host divergence

**Setup:** `envs/panel.yaml` (mmseqs2, samtools, seqkit),
`workflow/scripts/build_panel.py` (`build_panel` rule) fully implements
Section 8:
1. Adjacent windows of the same region class from `panel_regions.bed`
   (7.3) merged into contiguous blocks (only panel-relevant classes:
   strict_core/soft_core/shell/private_accessory/accessory_chromosome/
   starship_like - repeat_ambiguous/unclassified/subtelomeric_dynamic
   deliberately excluded, as too uncertain for panel anchors).
2. Sequences of all blocks extracted across all 5 genomes via
   `samtools faidx -r` (4,213 blocks).
3. **A single mmseqs2 easy-cluster run (98% identity, 90% mutual
   coverage, Section 8.2) on ALL blocks together** implements both
   halves of the dedup rule simultaneously: near-identical copies of ONE
   homologous block across multiple genomes collapse into a single
   cluster/representative; structurally different variants (even among
   core/soft-core/shell) remain automatically separated - no special
   handling per region type needed.
4. Panel manifest (`results/panel/panel_contig_manifest.tsv`, exactly
   the document's columns) + final panel FASTA
   (`results/panel/Mo_multiref_panel_v1.fa` + `.fai`) with
   `>PANEL######|type=...|cluster=...|rep=...|source=...` headers
   (Section 8.3), cluster prefixes `CORE_`/`ACC_`/`STAR_` depending on
   region type (Section 8.4 example).

**Bug (duplicated header prefix):** The first run failed with a
`KeyError` - `block_id()` already returns
`"{genome_id}__{contig}:{start}-{end}"` (since `contig` itself already
carries the genome prefix), but the lookup key was incorrectly built as
`(genome_id, block_id(b))` - a duplicated prefix. Fix: lookup directly
via the full `block_id()` string, which matches the FASTA/mmseqs2
sequence name exactly.

**Resources:** mmseqs2's prefiltering index build briefly needed ~9.5 GB
of the 10 GB WSL RAM (223 MB free) - close to the OOM limit but
survived; memory demand then dropped to ~5 GB for the actual search/
alignment phase. Total runtime about 5 minutes for 7,582 extracted
sequences (considerably more than originally expected, see the next
point).

**Important methodological finding:** The dedup rate is considerably
lower than the document's example suggests (there: ONE representative
covers up to 14 genomes). Result here: **3,677 panel regions from 4,213
pre-blocks** - only ~13% collapse. Of 1,783 `strict_core` clusters
(gene-based, present in ALL 5 genomes), **1,590 (89%) contain only ONE
genome** in the cluster - the corresponding blocks of the other 4
genomes are, despite orthogroup homology, too divergent (>2% sequence
divergence at the 10-kb window level) to cluster at 98% identity.

**Explanation:** This is a direct, expected consequence of the panel's
own design decision (5 maximally divergent host representatives instead
of close relatives, see the "panel reduced to 5 host representatives"
entry). At ~1 SNP per 160 bp between host lineages (from the SyRI
figures: 237k-272k SNPs across ~43 Mb), the expected identity even in
genuine ortholog blocks between two host lineages lies in the 98-99%
range only for short, low-variability stretches - a 10-kb window often
just misses the 98% threshold.

**Assessment - not a bug, but a trade-off:** The resulting panel is
larger than the document's example, but more informative as a result:
it retains the actual lineage-specific sequence variants of "core"
regions, instead of replacing them with a single (e.g., Oryza-only)
representative that would be less suitable for long-read mapping from
other host lineages (exactly the problem Section 7.2 warns about
regarding a single anchor). The 98%/90% threshold remains unchanged at
the document's value (`thresholds.yaml: panel.dedup_identity/
dedup_coverage`) - lowering it would change the panel's biological
meaning (more compression, but loss of lineage-specific core variants)
and is deliberately NOT done here, but left as a configurable,
documented parameter.

**Panel composition (3,677 regions, 184.5 MB FASTA):**

| Region type | Panel clusters |
|---|---|
| `strict_core` | 1,783 |
| `soft_core` | 1,610 |
| `shell` | 245 |
| `starship_like` | 30 |
| `accessory_chromosome` | 5 |
| `private_accessory` | 4 |
