# Analysis Plan: Starfish, Starship, and Mini-Chromosome Comparison

**Project:** PhD thesis – Barragan / *Magnaporthe oryzae*  
**Target audience:** Entry-level, without extensive bioinformatics background  
**Date:** August 31, 2026

---

## 1. Overview: What is being built?

You are not developing a single analysis, but a **reproducible pipeline**. A pipeline is a fixed, documented sequence of steps: every input file, every parameter, and every generated result file is defined. This means the analysis can later be applied to 1, 40, or 417 isolates without having to manually copy commands.

The pipeline answers three distinct but interconnected questions:

| Analysis area | Central question | Typical evidence |
|---|---|---|
| Starfish / Starships | Are there YR-associated large mobile elements? | YR HMM hits, Captain candidates, gene context, cargo genes, boundaries |
| Mini-chromosomes | Which sequences are likely accessory mini-chromosomes? | Size, repeat content, core-gene scarcity, coverage, synteny, homology |
| Comparison/HGT | Which elements are unusually similar or shared between isolates? | PAV, alignment coverage, sequence identity, core-vs-element contrast, phylogenetic discordance |

An important principle:

> A YR hit is not automatically a Starship. A short contig is not automatically a mini-chromosome. High sequence similarity is not automatically horizontal transfer.

For robust conclusions, you need multiple independent lines of evidence.

---

## 2. Key terms

### Assembly

An **assembly** is the reconstructed genome sequence of an isolate. It consists of one or more sequences, called contigs or chromosomes.

### Contig

A **contig** is a contiguously assembled DNA sequence. In a very good long-read assembly, a contig can correspond to a complete chromosome. In a fragmented assembly, a chromosome can be distributed across multiple contigs.

### GFF

A **GFF file** (*General Feature Format*) describes where genes and other features are located on the DNA. It contains, among other things:

- Sequence name or chromosome/contig
- Start and end position of a gene
- Strand direction
- Feature type, e.g. `gene`, `mRNA`, `CDS`
- IDs and further attributes

The contig ID in the first GFF column must exactly match the header of the corresponding FASTA.

### YR and Captain gene

YR stands for **tyrosine recombinase**. Tyrosine recombinases can mediate DNA recombination. In Starships, a YR gene is often a Captain gene associated with the mobility of the element. Not every YR gene is necessarily part of a Starship.

### Starship

A **Starship** is a large mobile genetic element that carries a Captain gene and, variably, additional genes often referred to as cargo. To convincingly annotate a Starship, more than the detection of a YR domain is needed: gene environment, plausible boundaries, repeats, and a plausible element structure are important.

### Mini-chromosome

A **mini-chromosome** or accessory chromosome is a chromosomal unit that does not belong to the conserved core genome. Such sequences can differ substantially between isolates and can contain genes involved in adaptation, pathogenicity, or horizontal transfer.

---

## 3. Why Snakemake?

**Snakemake** is a workflow engine for bioinformatics. You define rules instead of executing a long list of manual commands.

Example concept:

```text
normalized FASTA
        ↓
Starfish YR annotation
        ↓
consolidated annotation
        ↓
Starship candidates
        ↓
PAV and similarity matrices
```

Snakemake recognizes which intermediate products are missing and only reruns the necessary steps. It can process multiple isolates in parallel and stores logs per rule.

### Advantages

- Identical procedure for all isolates
- No manual copying of commands
- Resumable after errors
- Parallel execution
- Clear dependencies
- Good foundation for supplement, thesis, and manuscript
- Exact documentation of the files and parameters used

### Why not one large Python script?

Python is useful for custom calculations and tables. However, a single large script quickly becomes hard to maintain: it must itself check which files exist, which steps need to rerun, and how parallel jobs are coordinated. Snakemake handles precisely this orchestration. Python is used as a tool within individual pipeline steps.

---

## 4. Recommended project structure

Set up the analysis as clearly separated directories:

```text
barragan/
├── config/
│   ├── samples.tsv
│   ├── paths.yaml
│   ├── parameters.yaml
│   └── references.yaml
├── workflow/
│   ├── Snakefile
│   ├── rules/
│   │   ├── qc.smk
│   │   ├── headers.smk
│   │   ├── starfish.smk
│   │   ├── annotation.smk
│   │   ├── minichromosomes.smk
│   │   ├── comparisons.smk
│   │   └── reporting.smk
│   ├── scripts/
│   │   ├── normalize_fasta_headers.py
│   │   ├── normalize_gff_seqids.py
│   │   ├── parse_starfish_results.py
│   │   ├── calculate_contig_metrics.py
│   │   ├── define_mchr_candidates.py
│   │   ├── build_pav_matrix.py
│   │   └── score_hgt_candidates.py
│   └── envs/
│       ├── starfish.yaml
│       ├── alignment.yaml
│       └── python.yaml
├── input/
│   ├── assemblies/
│   ├── annotations/
│   └── metadata/
├── results/
│   ├── qc/
│   ├── normalized/
│   ├── starfish/
│   ├── minichromosomes/
│   ├── comparisons/
│   └── figures/
├── logs/
├── docs/
│   ├── workflow.md
│   ├── decisions.md
│   └── provenance.md
└── README.md
```

### What belongs where?

| Location | Content |
|---|---|
| `input/` | Raw data or controlled inputs; not to be silently modified |
| `config/` | Tables and parameters that control the run |
| `workflow/` | Snakemake rules, small analysis scripts, and software environments |
| `results/` | Final and intermediate results produced by the pipeline |
| `logs/` | Complete standard output and error messages per rule/isolate |
| `docs/` | Decisions, method descriptions, and data provenance |

---

## 5. Build metadata first

Before running many analyses, each isolate needs a unique, stable identifier. Create `config/samples.tsv` with at least:

```tsv
isolate_id	assembly_accession	assembly_fasta	gff	host	lineage	country	year	sequencing	include
ASM434696v1	GCA_004346965.1	/path/to/GCA_004346965.1_ASM434696v1_genomic.fna	/path/to/GCA_004346965.1_ASM434696v1_genomic.gff	unknown	unknown	unknown	unknown	NCBI	TRUE
```

### Meaning of the key columns

| Column | Purpose |
|---|---|
| `isolate_id` | Short, unique ID for tables, figures, and file names |
| `assembly_accession` | NCBI accession identifier, e.g. `GCA_004346965.1` |
| `assembly_fasta` | Full path to the original assembly |
| `gff` | Full path to the gene annotation; empty if not available |
| `host`, `lineage`, `country`, `year` | Biological metadata for later interpretation |
| `sequencing` | e.g. Illumina, PacBio, Nanopore, or unknown |
| `include` | Controls whether an isolate is used in a given run |

**Rule:** An ID must never be changed "silently" after the fact. If an ID is corrected, document the reason in `docs/provenance.md`.

---

## 6. Phase A: Checking and standardizing inputs

### A1. Capture assembly quality

First, a quality overview is generated per assembly. At minimum, capture:

- Number of contigs
- Total size
- Longest contig
- N50
- Number of `N` bases
- GC content
- Available annotation
- Sequencing technology, if known

**Why?** Differences in contig count and fragmentation can mimic later differences in mini-chromosomes. A missing mChr may be biologically absent — or simply not assembled.

### A2. Make FASTA headers unique

Starfish expects, or benefits from, headers in the scheme:

```text
>genomeID_contigID
```

Example:

```text
Original:      >CP034210.1
Normalized:    >GCA_004346965.1_CP034210.1
```

This is important for the multi-isolate analysis, because contig IDs would otherwise not be unique across isolates.

### A3. Adjust the GFF in sync

If the FASTA headers are changed, the first column of the GFF must be changed synchronously.

Example:

```text
Original GFF column 1:      CP034210.1
Normalized GFF column 1:    GCA_004346965.1_CP034210.1
```

**Quality control:** After the adjustment, all seqIDs from the GFF must be present in the FASTA. Otherwise, the FASTA and GFF must not be used together.

### A4. Leave original data unchanged

The NCBI files remain unchanged in a raw-data folder. The pipeline writes normalized copies to, for example:

```text
results/normalized/{isolate_id}.fna
results/normalized/{isolate_id}.gff
```

This makes it possible to trace what was changed at any time.

---

## 7. Phase B: Starfish workflow

### B1. Finding YR candidates

The first Starfish step searches for YR genes de novo. For this you use:

```text
YRsuperfamRefs.faa
YRsuperfams.p1-512.hmm
```

The pattern is:

```bash
starfish annotate \
  --assembly <assemblies.tsv> \
  --profile <YRsuperfams.p1-512.hmm> \
  --proteins <YRsuperfamRefs.faa> \
  --prefix <isolate>_YR \
  --idtag YR \
  --outdir <output_dir> \
  --tempdir <temp_dir> \
  --threads <n>
```

### What happens internally?

1. MetaEuk searches the DNA for possible protein-coding genes similar to the reference YR proteins.
2. `hmmsearch` checks the predicted proteins against the YR HMM.
3. Hits that pass the HMM threshold are output as filtered candidates.

### B2. Naming result statuses correctly

Use clear status names:

| Status | Definition |
|---|---|
| `YR_candidate` | HMM-validated YR gene from `starfish annotate` |
| `contextual_Starship_candidate` | YR candidate with a plausible genomic context |
| `bounded_Starship_candidate` | Candidate with additionally plausible boundaries/flanking features |
| `high_confidence_Starship` | Predefined criteria for structure and evidence fulfilled |

This terminology prevents a preliminary hit from being over-interpreted.

### B3. Incorporating existing gene annotation

If NCBI GFF files are available and compatible with the normalized FASTA, pass them via `-g`:

```text
genomeID<TAB>path-to-normalized-GFF
```

Why? The de novo YR search primarily finds candidates for YR genes. For cargo genes and neighborhood analyses, you need the most complete gene models possible across the whole genome.

### B4. Subsequent Starfish steps

The local Starfish installation should be documented with `starfish --help`. Typical steps in the Starfish workflow are:

1. Annotation or consolidation of the existing and new gene models.
2. Determination of the gene context around YR candidates.
3. Search for insertion sites and/or boundaries.
4. Analysis of flanking sequences and possible repeats.
5. Summary of the candidates and possible cargo genes.

**Important:** The exact subcommands, parameters, and file formats are documented from the locally installed version before implementation. Not all Starfish versions have the same names or options.

### B5. Starfish output per isolate

Recommended final table:

```text
results/starfish/starship_candidates.tsv
```

Possible columns:

```tsv
candidate_id	isolate_id	contig_id	start	end	strand	yr_gene_id	yr_hmm_evalue	yr_hmm_score	context_class	boundary_class	cargo_gene_count	cargo_annotations	confidence	comments
```

---

## 8. Phase C: Finding mini-chromosomes

Mini-chromosomes must be investigated independently of Starfish. Starfish is not a general-purpose tool for reliable mChr detection.

### C1. Calculate metrics per contig

Create a table per contig, e.g.:

```text
results/minichromosomes/contig_metrics.tsv
```

Recommended columns:

```tsv
isolate_id	contig_id	length_bp	gc_fraction	gene_count	core_gene_count	repeat_fraction	TE_fraction	secreted_protein_count	effector_candidate_count	median_depth	depth_ratio_to_core	mchr_score	classification
```

Not all columns are available from the start. Begin with length, GC content, and gene count; expand step by step afterward.

### C2. Evidence for mChr candidates

A contig is not labeled as an mChr based on a single threshold. Use a combination of features:

- small size relative to the core chromosomes
- low density of conserved core genes
- high repeat or transposon content
- high density of effectors or secreted proteins, if biologically plausible
- conspicuous coverage or copy number, if raw reads are available
- low synteny to core chromosomes
- strong homology to known accessory regions or mChrs

### C3. Working classes

Use conservative classes instead of a premature yes/no decision:

| Class | Meaning |
|---|---|
| `core_like` | apparently part of a conserved core chromosome |
| `accessory_candidate` | shows one or more accessory properties |
| `mChr_candidate` | multiple lines of evidence support a mini-chromosome |
| `uncertain` | assembly fragmentation or evidence does not allow a reliable classification |

### C4. Dealing with fragmented assemblies

An mChr can be distributed across multiple contigs in a poor assembly. Therefore, the analysis must not exclusively ask: "Which single contig is an mChr?"

Homologous sequence blocks should also be considered. Several small contigs of the same isolate can jointly represent parts of the same mChr if they align to the same reference mChr or the same mChr cluster from different isolates.

---

## 9. Phase D: Comparison between isolates

### D1. Determining the core genome as background

Before HGT is interpreted, a benchmark for comparison is needed: How similar are the isolates in the conserved genome?

Possible approaches:

- Determine and align single-copy orthologs
- Identify conserved core regions
- Calculate pairwise distances from core SNPs or core alignments

Output example:

```text
results/comparisons/core_distance.tsv
```

### D2. Comparing mChr and Starship sequences

First, a quick pre-filter can group similar sequences. This is followed by a precise comparison via sequence alignment.

Recommended principles:

- consider both alignment directions
- store query and target coverage
- store the percentage identity per pair
- document oriented and collinear blocks
- do not rely only on the best short hit

Output example:

```text
results/comparisons/pairwise_element_alignments.tsv
```

Possible columns:

```tsv
query_isolate	query_element	target_isolate	target_element	query_coverage	target_coverage	aligned_bp	identity	orientation	synteny_class
```

### D3. Presence/absence variation (PAV)

PAV means: Is an element or a homologous block present, absent, or unclear in an isolate?

Example of a PAV matrix:

| Element cluster | Isolate_A | Isolate_B | Isolate_C |
|---|---:|---:|---:|
| mChr_cluster_001 | present | present | absent |
| Starship_cluster_014 | present | uncertain | absent |

In the case of fragmentation, `uncertain` is preferable to a false `absent`.

### D4. Core-vs-element contrast

A central HGT logic states:

> If two isolates are relatively different in the core genome but share an mChr or Starship that is nearly identical and extensively conserved, this is an indication of a more recent exchange or transfer of the element.

Example illustration:

```text
Isolate A vs. B
Core genome identity:                  97.5 %
Identity of an mChr cluster:           99.95 %
Reciprocal mChr coverage:              92 %
```

This example is not a fixed decision rule. The thresholds must be calibrated based on your assemblies, your isolate population, and, where possible, known positive and negative comparison cases.

### D5. Phylogenetic discordance

An additional strong signal is when the ancestry relationship of an element does not match the core-genome relatedness.

Example:

- In the core genome, isolates A and B do not group closely together.
- Their mChr sequences, however, group closely and show high, extensive sequence similarity.

This can be consistent with a transfer. However, alternatives such as conserved selection, contamination, misassembly, or incomplete sampling must be checked.

---

## 10. Minimal first pipeline milestone

Do not build a complete large workflow first. Reach these verifiable milestones in sequence.

### Milestone 1: Fully prepare one isolate

- Create `samples.tsv`
- Locate the original FASTA and, if applicable, the original GFF
- Normalize FASTA headers
- Normalize GFF seqIDs in sync
- Verify that all GFF seqIDs exist in the FASTA

**Success criterion:** The normalized FASTA and GFF match and are documented.

### Milestone 2: Run Starfish stably on one isolate

- Repeat the YR annotation with the normalized contig IDs
- Check the log for header warnings and tool errors
- Inspect the filtered YR GFF
- Convert the result into a tabular candidate list

**Success criterion:** A reproducible Starfish output with unique sequence and gene IDs.

### Milestone 3: Analyze three isolates as a pilot

Choose, if possible:

- one isolate with known or expected mChr evidence
- one closely related isolate
- one more distantly related isolate

Compare contig metrics, Starfish candidates, and sequence homologies.

**Success criterion:** The pipeline runs multiple times without manually changing commands; differences and limitations of the data become visible.

### Milestone 4: Scale to all isolates

Only once the pilot is plausible:

- activate all isolates in `samples.tsv`
- set available threads and memory realistically
- generate logs and checks per step
- merge results into overall tables

---

## 11. Quality controls and sources of error

### Missing or fragmented mChrs

If an assembly consists only of short contigs, an existing mChr can be fragmented or not fully assembled at all. "Not found" is then not necessarily "biologically absent."

### Different sequencing quality

Long-read assemblies are usually better suited than heavily fragmented short-read assemblies for detecting complete mChrs, repeats, and element boundaries. Assembly quality must therefore be carried along in every comparison table.

### Contamination and misassignment

Unusually similar elements between distantly related isolates can indicate HGT, but can also arise from contamination, sample mix-up, or misassigned contigs. For strong candidates, check read mapping, coverage, gene taxonomy, and assembly provenance.

### Over-interpretation of a YR hit

A YR HMM hit is the beginning, not the end, of the analysis. Maintain the distinction between YR candidate, Starship candidate, and high-confidence element in all tables and figures.

### Changing thresholds after viewing results

Note criteria and parameters early in `docs/decisions.md`. If changes later become necessary, version them with justification and repeat the affected analyses.

---

## 12. Documentation during the work

Create and maintain at least three documents.

### `README.md`

Contains:

- Goal of the pipeline
- Brief installation instructions
- Example call
- Folder overview
- Expected main outputs

### `docs/workflow.md`

Contains:

- graphical or textual sequence of the rules
- Input and output per rule
- Tools used
- Expected runtime/resource order of magnitude

### `docs/decisions.md`

Contains methodological decisions, e.g.:

```text
2026-08-31
Decision: FASTA headers are normalized as <assembly_accession>_<original_contig_id>.
Rationale: Unique contig IDs for multi-isolate analyses and compatibility with Starfish.
Consequence: Associated GFF seqIDs are adjusted synchronously.
```

Also document:

- the Starfish version used
- database files and HMM thresholds
- MetaEuk parameters
- criteria for mChr classes
- criteria for Starship confidence
- criteria for shared elements and HGT prioritization

---

## 13. Concrete next steps

The most sensible order is:

1. **Inventory the raw data.** For each isolate, capture FASTA, GFF, possible raw reads, and metadata in `samples.tsv`.
2. **Program the normalization.** Create a small Python script that safely renames FASTA headers and writes a mapping table `old_id → new_id`.
3. **Program the GFF synchronization.** A second script that uses the mapping table and replaces the seqIDs in column 1.
4. **Build in validation.** A rule that checks that GFF seqIDs and FASTA seqIDs match.
5. **Repeat the Starfish test.** Test `GCA_004346965.1` with the normalized FASTA to eliminate the previous header warnings.
6. **Document the subsequent Starfish workflow locally.** Save `starfish --help` and the help for all subcommands used; only then implement the exact steps as Snakemake rules.
7. **Generate contig metrics.** Calculate length, GC content, and gene count for all contigs; this is the starting point of the mChr analysis.
8. **Plan the pilot comparison.** Select three representative isolates and comparatively align their mChr candidates as well as YR/Starship candidates.
9. **Only then scale up.** After a successful pilot run, process the entire isolate collection.

---

## 14. Practical guideline

Work iteratively and keep every step small:

```text
1 isolate → 3 isolates → entire dataset
1 rule    → test output → documentation → next rule
```

If a result is surprising, first check the inputs, headers, assembly quality, and mapping coverage. A cleanly documented negative or uncertain result is scientifically more valuable than a quick but unverifiable claim.
