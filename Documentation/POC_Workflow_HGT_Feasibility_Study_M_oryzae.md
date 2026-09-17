# POC Workflow: Feasibility Study HGT in *M. oryzae*

Guideline for technical implementation in VS Code · Reference: dissertation proposal "Population Genomics of Horizontal Gene Transfer in Magnaporthe oryzae", WP1.2–WP1.4, Q1–Q3 · As of 2026\-09\-01

## Purpose of This Document

This document is a **working basis for VS Code**, not a final result. It defines a small-scale proof\-of\-concept (POC) intended to show, before application submission/funding approval, that the methodology described in the proposal (WP1.2–WP1.4) works technically and that the planned analysis (Q3, GLMM) is statistically sound. It covers the four requested points:

1. Functionality of the multi-reference panel
2. Number of Starships and other HGT candidates in the POC
3. Clustering of core vs. candidate regions by host type
4. Power analysis for the planned GLMM (Q3)

**Important assumption:** No real sequence data are available in this session (only the proposal PDF in the project). Steps 1–3 are therefore formulated as an **executable protocol** with concrete tools/commands/parameters, ready to apply as soon as reference genomes and an isolate subset (e.g. from the 417\-isolate dataset, Barragan et al. 2024, Suppl. Table S5) are available. Step 4 (power analysis) is data-independent and can already be simulated now with assumed effect sizes — these assumptions are explicitly marked and must be updated as soon as initial pilot data from WP1.3 are available.

## Recommended Repo Structure

```text
poc-hgt-magnaporthe/
├── data/
│   ├── references/            # Multi-reference panel (FASTA + GFF), see Section 1
│   └── isolates_poc/          # Illumina/Nanopore subset for the POC
├── envs/                      # conda/mamba environment.yml per tool
├── workflow/
│   ├── Snakefile              # or nextflow main.nf
│   ├── rules/                 # mapping.smk, pav.smk, starfish.smk, cluster.smk
│   └── scripts/
├── results/
│   ├── mapping/
│   ├── pav_calls/
│   ├── starfish/
│   └── clustering/
├── power_analysis/
│   └── glmm_power_sim.R       # see Section 4
└── README.md                  # ← this document (as export)
```

A workflow manager (Snakemake or Nextflow) is recommended, since WP1.3 explicitly calls for a "modular pipeline" with calibratable thresholds — reproducible and later directly scalable to the full 417\-isolate dataset as well as the CAU HPC cluster.

* * *

## 1\. Functionality of the Multi-Reference Panel

### Goal

Test whether the mapping procedure described in the proposal (short reads against a multi-reference panel, coverage-/PAV-based candidate detection) works technically and whether the proposed thresholds (80–90% breadth at ≥5× depth = "present", \<10–20% breadth = "absent") are calibratable.

### Reference Selection (POC Scope: 4–6 Genomes)

Criterion according to the proposal: coverage of several host lineages plus, where possible, high-quality, published long-read assemblies.

| Lineage / Host | Reference | Source |
| --- | --- | --- |
| Rice (*Oryza sativa*) | 70\-15 (near\-complete) | Cheng et al. 2025 |
| Wheat (*Triticum aestivum*) | Wheat\-blast reference isolate | Yoshida et al. 2016 |
| Finger millet (*Eleusine coracana*) | Eleusine lineage isolate | Chiapello et al. 2015 |
| Wild-grass lineages | ≥1 isolate from Barragan et al. 2022/2024 | CAU / bioRxiv supplement |
| *M. grisea* (sister taxon, negative control) | Reference assembly | Chiapello et al. 2015 |

### Pipeline Steps

```bash
# 1. QC of the short reads
fastp -i isolate_R1.fq.gz -I isolate_R2.fq.gz \
      -o qc_R1.fq.gz -O qc_R2.fq.gz --json qc.json

# 2. Index panel (one index per reference, or combined panel)
bwa-mem2 index references/panel_combined.fa
# Alternative: minimap2 -d panel.mmi references/panel_combined.fa

# 3. Mapping
bwa-mem2 mem -t 8 references/panel_combined.fa qc_R1.fq.gz qc_R2.fq.gz \
  | samtools sort -@4 -o mapping/isolate.sorted.bam
samtools index mapping/isolate.sorted.bam

# 4. Coverage / breadth per candidate region
mosdepth --by candidate_regions.bed -t 4 pav_calls/isolate mapping/isolate.sorted.bam

# 5. PAV call (script, thresholds parameterized)
python scripts/pav_call.py \
  --regions candidate_regions.bed \
  --mosdepth pav_calls/isolate.regions.bed.gz \
  --present-breadth 0.85 --present-depth 5 \
  --absent-breadth 0.15 \
  --out pav_calls/isolate.pav.tsv
```

### Calibration & Success Criteria

- Calibration of the thresholds against **known mini-chromosome regions** (mChrA, Barragan et al. 2024) and replicate samples.
- Goal: **sensitivity ≥ 90%, specificity ≥ 95%**, final thresholds via F1-score maximization (grid search over breadth/depth combinations).
- Sensitivity analysis: vary thresholds in 5% steps, validate PAV calls against long-read-confirmed regions.

### POC Feasibility Questions (to Be Answered)

- [ ] Do reads from divergent lineages (e.g. a wild-grass isolate against the rice reference) still map reliably enough for a coverage signal?
- [ ] How strongly does reference bias (which reference in the panel "wins") affect the PAV calls? → Test: same isolate against a single reference vs. combined panel.
- [ ] Does the workflow scale computationally to the full dataset (runtime/memory per isolate × 417)?

* * *

## 2\. Starship and HGT Candidate Counting in the POC

### Goal

For the POC subset (see above, 4–6 references + 10–20 isolates), determine **how many Starships and other HGT candidate regions** are found using the two complementary methods described in the proposal.

### Two Complementary Detection Approaches (per Proposal, Method a+b)

**a) Coverage-based (multi-reference panel, see Section 1):** PAV regions from `pav_call.py`, filtered for candidates with discordant clustering (see Section 3).

**b) Structure-based (Starfish / Stargraph):**

```bash
# Starfish: detection of Starship elements based on Captain tyrosine recombinase
# + characteristic flanking signatures (Gluck-Thaler & Vogan 2024)
starfish annotate -T 8 -x isolate_id \
  -a references/isolate_assembly.fa \
  -g references/isolate_assembly.gff3 \
  -o results/starfish/isolate_id

starfish insert -x isolate_id \
  -a results/starfish/isolate_id.starships.bed \
  -o results/starfish/isolate_id.inserts

# Stargraph: pan-genomic graph approach for consistency checking across isolates
stargraph build --assemblies references/*.fa --starships results/starfish/*.bed \
  -o results/starfish/pangenome_graph
```

### Cross-Validation (per Proposal)

- Coverage-based PAV candidates → screen for Starship hallmarks (Captain protein, terminal repeats, target site duplications).
- Starfish candidates → check for coverage discordance across lineages.
- Classification of each candidate region into a **vector class**: `mini-chromosome`, `Starship`, `both`, `unclear`.

### Expected POC Output

A candidate table as the central intermediate result:

| Field | Description |
| --- | --- |
| `region_id` | Unique ID of the candidate region |
| `vector_class` | mChr / Starship / both / unclear |
| `size_kb` | Size (Starships: 25–700 kb per definition) |
| `n_isolates_present` | Number of POC isolates with a "present" call |
| `host_lineages_present` | Affected host lineages (for Section 3) |
| `captain_detected` | yes/no (Starships only) |
| `evidence_tier` | high-confidence / exploratory (per risk assessment Q1) |

**Target for the POC:** Since actual numbers are only available after the pipeline has run, the POC success metric is not "N Starships" but rather: *Is at least 1 known, published Starship region (O'Donnell et al. 2025) correctly recovered, and does the pipeline yield a non-trivial but manageable number of candidates (rough expectation in the low double-digit range for ~10–20 isolates, no fixed target)?* This serves as a sensitivity sanity check, not as a preprint figure.

* * *

## 3\. Clustering of Core vs. Candidate Regions by Host Type

### Goal

Reproduce the logic shown in **Figure 2** of the proposal: core genome markers should cluster isolates by host-associated lineage, whereas candidate HGT regions should show **discordant** (cross-lineage) clustering — the central evidence criterion for HGT (Q1).

### Methodology

1. **Core genome matrix:** PAV or SNP matrix from core genome-wide markers (outside the candidate regions).
2. **Candidate matrix:** PAV matrix using only the candidate regions identified in Section 2.
3. For both matrices: PCA (or MDS), colored by `host_lineage` (rice, wheat, Eleusine, wild grass).
4. **Quantify discordance:** compare clusters from PCA (k-means or hierarchical clustering, k = number of known lineages) with the known host labels via the **Adjusted Rand Index (ARI)**.
   - Expectation for the core matrix: **high ARI** (clusters ≈ host lineage).
   - Expectation for the candidate matrix: **lower ARI**, with clearly identifiable cross-lineage groups (= HGT signature).

```python
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

def cluster_and_score(pav_matrix: pd.DataFrame, host_labels: pd.Series, k: int):
    pca = PCA(n_components=min(10, pav_matrix.shape[1]))
    coords = pca.fit_transform(pav_matrix.values)
    clusters = KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(coords)
    ari = adjusted_rand_score(host_labels, clusters)
    return coords, clusters, ari

core_coords, core_clusters, core_ari = cluster_and_score(core_pav, host_labels, k=4)
cand_coords, cand_clusters, cand_ari = cluster_and_score(candidate_pav, host_labels, k=4)
print(f"Core ARI: {core_ari:.2f}  |  Candidate ARI: {cand_ari:.2f}")
```

### Recommendation: Method Validation on Synthetic Data (Already Possible Now)

Before real data are available, the analysis logic itself can be validated on a synthetic toy dataset (known "ground truth" HGT isolates artificially inserted, expectation: the pipeline must correctly detect the injected discordance). A script for this (`scripts/validate_clustering_synthetic.py`) is recommended as a standalone test that checks:

- Core matrix without injected HGT events → ARI close to 1.
- Candidate matrix with artificially inserted cross-lineage sharing (e.g. 15% of isolates share a "foreign" candidate region) → ARI significantly lower, affected isolates visibly grouped in the candidate PCA.

This serves as a **unit test of the analysis method**, not as a biological result.

### Interpretation Criteria

- [ ] Core ARI clearly \> candidate ARI (qualitative, no fixed threshold specified in the proposal)
- [ ] Visual discordance reproduces the pattern from Fig. 2 (wild-grass/rice isolates cluster together in candidate regions)
- [ ] Result consistent across multiple k values / clustering methods (robustness check)

* * *

## 4\. Power Analysis for the GLMM (Q3)

### Goal

Estimate whether a POC/pilot sample (or the later full field study) has sufficient power to detect the effects specified in the proposal: Poisson or negative-binomial GLMM with fixed effects (donor presence, recipient presence, vector type, donor×vector and recipient×vector interaction), random intercepts for field and year, offset = log(number of isolates sampled per field-year).

### Recommended Approach: Simulation-Based Power Analysis (Monte Carlo)

**Recommended tooling:** R with `glmmTMB` (or `lme4::glmer.nb`) + `simr` — standard tools for exactly this model in ecology/evolutionary biology (cf. Bolker et al. 2009, cited in the proposal). *Note: R is not installed in this cloud session — this section provides a code skeleton directly runnable in VS Code; as a fallback without R, `statsmodels.genmod.bayes_mixed_glm.PoissonBayesMixedGLM` in Python is suitable, with lower accuracy for small cluster numbers.*

### POC Sample Scenarios (Assumptions — to Be Updated with Pilot Data from WP1.3)

| Scenario | Fields | Years | Isolates / Field-Year | Total Isolates | Assumed HGT Baseline Rate | Assumed Donor×Vector Effect (Rate Ratio) |
| --- | --- | --- | --- | --- | --- | --- |
| Conservative (minimal POC) | 5 | 2 | 10 | 100 | 5% | 1\.5× |
| Moderate (realistic POC) | 10 | 2 | 20 | 400 | 8% | 1\.8× |
| Optimistic (extended POC) | 15 | 3 | 30 | 1\.350 | 12% | 2\.0× |
| Full study (reference) | \>20 | ≥3 | variable | \~1,000\+ | unknown | — |

These values are **placeholders**, not published effect sizes — HGT event rates are, per the proposal itself (risk assessment Q3), still unknown and potentially rare/unevenly distributed.

### Code Skeleton (`power_analysis/glmm_power_sim.R`)

```r
library(glmmTMB)
library(simr)

set.seed(1)

simulate_scenario <- function(n_fields, n_years, n_iso_per_fy,
                               baseline_rate, effect_rr) {
  df <- expand.grid(field = factor(1:n_fields), year = factor(1:n_years))
  df$log_n_isolates <- log(n_iso_per_fy)
  df$donor     <- rbinom(nrow(df), 1, 0.5)
  df$recipient <- rbinom(nrow(df), 1, 0.5)
  df$vector    <- factor(sample(c("mChr", "Starship", "both"),
                                 nrow(df), replace = TRUE))

  beta <- c(intercept = log(baseline_rate),
            donor = 0.3, recipient = 0.2,
            vectorStarship = 0.2, vectorboth = 0.4,
            donor_vectorStarship = log(effect_rr))

  # linear predictor including random intercepts (field, year) + Poisson draw
  field_re <- rnorm(n_fields, 0, 0.4)[df$field]
  year_re  <- rnorm(n_years, 0, 0.3)[df$year]
  lp <- beta["intercept"] + beta["donor"] * df$donor +
        beta["recipient"] * df$recipient +
        ifelse(df$vector == "Starship", beta["vectorStarship"],
               ifelse(df$vector == "both", beta["vectorboth"], 0)) +
        ifelse(df$donor == 1 & df$vector == "Starship",
               beta["donor_vectorStarship"], 0) +
        field_re + year_re + df$log_n_isolates

  df$hgt_count <- rpois(nrow(df), exp(lp))
  df
}

power_for_scenario <- function(..., nsim = 200) {
  df <- simulate_scenario(...)
  fit <- glmmTMB(hgt_count ~ donor * vector + recipient * vector +
                   offset(log_n_isolates) + (1 | field) + (1 | year),
                 family = poisson, data = df)
  powerSim(fit, test = fixed("donor:vectorStarship"), nsim = nsim)
}

# Run through scenarios from the table above and log power (%)
scenarios <- list(
  konservativ  = list(n_fields = 5,  n_years = 2, n_iso_per_fy = 10, baseline_rate = 0.05, effect_rr = 1.5),
  moderat      = list(n_fields = 10, n_years = 2, n_iso_per_fy = 20, baseline_rate = 0.08, effect_rr = 1.8),
  optimistisch = list(n_fields = 15, n_years = 3, n_iso_per_fy = 30, baseline_rate = 0.12, effect_rr = 2.0)
)

results <- lapply(scenarios, function(s) do.call(power_for_scenario, s))
```

### Interpretation & Feasibility Decision

- **Rule of thumb:** With so few clusters (5–15 fields × 2–3 years), the power for **interaction effects** (donor×vector) is typically markedly lower than for main effects — this is the most critical test in the simulation.
- If even the "optimistic" POC scenario shows \< 80% power for the interaction effect: (a) simplify the model (interaction only exploratory, as foreseen as a fallback in the proposal: "stepwise, beginning with simple models"), or (b) size the sample for the full study accordingly.
- Once WP1.3 delivers initial real HGT candidate rates, replace `baseline_rate` and `effect_rr` with observed values and repeat the simulation — the present analysis is a **starting point**, not a final proof.

### Feasibility Checklist

- [ ] Power ≥ 80% for main effects (donor, recipient, vector type) in the moderate scenario?
- [ ] Power ≥ 80% for at least one interaction effect in at least one realistic scenario?
- [ ] Sensitivity of power to random-effect variance (field/year) checked?
- [ ] Comparison of Poisson vs. negative binomial (overdispersion) planned, in case real count data are overdispersed?

* * *

## Summary: Open Points for Implementation in VS Code

- [ ] Arrange access to reference genomes (Section 1) and isolate subset (CAU HPC cluster)
- [ ] Set up environments: `bwa-mem2`/`minimap2`, `samtools`, `mosdepth`, `starfish`, `stargraph`, Python (`pandas`, `scikit-learn`), R (`glmmTMB`, `simr`)
- [ ] Create the Snakemake/Nextflow skeleton according to the repo structure above
- [ ] Implement and run the synthetic clustering validation test (Section 3) — possible immediately, independent of real data
- [ ] Run the power simulation (Section 4) in R locally/on HPC and fill in the scenario table with results
- [ ] After the first pipeline run: add the candidate table (Section 2) and ARI values (Section 3) to this document, update the assumptions in Section 4
