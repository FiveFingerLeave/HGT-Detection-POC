# Reproducible workflow: multi-reference panel and long-read PAV analysis in *Magnaporthe oryzae*

## Goal

This workflow builds a versioned, biologically annotated multi-reference panel from **14 complete genome assemblies**. It makes it possible to map long reads from test isolates against core, accessory, mini-chromosomal, and Starship-like genome regions, to determine presence/absence variation (PAV), and to assign each candidate region to one or more reference isolates.

Subsequently, **10 of 42 chromosome-level long-read isolates** from different host types are first analyzed as a pilot cohort. In parallel, a rarefaction analysis is carried out to test whether the reference panel already saturates the accumulation of additional candidate regions.

---

## 1. Research questions

1. Which sequence blocks in the 14-genome set are strict core, soft core, shell, or accessory?
2. Which complete or partial accessory chromosomes and mini-chromosomes occur in the panel?
3. Which regions are Starship-like, i.e. potentially mobile, cargo-bearing regions with DUF3435/Captain evidence?
4. Which of these regions are present, absent, or not reliably assessable due to repeats/multi-mapping in the ten long-read test isolates?
5. Which reference isolate(s) in the panel share a candidate region detected in a test isolate?
6. Does adding further reference genomes still yield many new candidate regions, or is the panel approaching saturation?

---

## 2. Basic principle

The workflow deliberately separates two reference products:

1. **Complete reference catalog:** All 14 genome assemblies remain complete, separate, and unchanged. It serves for reconstructing synteny, provenance, homology, and variants.
2. **Analytical multi-reference panel:** A non-redundant set of core representatives plus all validated alternative accessory regions, mini-chromosomes, and Starship-like regions. This panel is used for long-read mapping and quantitative PAV analysis.

Simply concatenating all 14 genomes is possible as an exploratory raw panel, but it is not ideal for quantitative PAV calls: similar core regions cause multi-mapping, spread coverage across multiple references, and can thus produce false absence calls.

---

## 3. Inputs

| Data type | Scope | Requirement |
|---|---:|---|
| Complete reference assemblies | 14 genomes | chromosome-/contig-level; including mini-chromosomes where possible |
| Gene annotation | 14 GFF3 files and proteomes or standardized reannotation | Uniform IDs and consistent pipeline |
| Repeat/TE annotation | 14 genomes | RepeatModeler2 + RepeatMasker or EDTA |
| Long reads | 42 isolates available; 10 pilot isolates | ONT or PacBio HiFi; metadata required |
| Starship references | DUF3435 HMM or known Captain proteins | For candidate search and annotation |
| Metadata | References and test isolates | Host, host group, origin, year, lineage, sequencing type, read N50, coverage |

---

## 4. Project structure

```text
magnaporthe_multiref_pav/
├── config/
│   ├── config.yaml
│   ├── references.tsv
│   ├── samples.tsv
│   └── thresholds.yaml
├── data/
│   ├── references/
│   ├── annotations/
│   ├── longreads/
│   └── resources/
├── envs/
│   ├── core.yaml
│   ├── annotation.yaml
│   ├── wga.yaml
│   ├── longreads.yaml
│   └── reporting.yaml
├── workflow/
│   ├── Snakefile
│   ├── rules/
│   │   ├── qc.smk
│   │   ├── annotation.smk
│   │   ├── repeats.smk
│   │   ├── starships.smk
│   │   ├── wga.smk
│   │   ├── panel.smk
│   │   ├── mapping.smk
│   │   ├── pav.smk
│   │   ├── rarefaction.smk
│   │   └── report.smk
│   └── scripts/
│       ├── rename_fasta_headers.py
│       ├── classify_regions.py
│       ├── build_panel.py
│       ├── call_pav.py
│       ├── annotate_pav.py
│       └── rarefaction.py
├── results/
├── logs/
├── README.md
└── LICENSE
```

Snakemake is suitable as a workflow engine because it supports rules, reproducible software environments, configuration files, parallel execution, logging, and the resumption of incomplete runs.

---

## 5. Metadata and unique IDs

### 5.1 Reference metadata

File: `config/references.tsv`

```tsv
genome_id	host	host_group	lineage	country	year	assembly_fasta	annotation_gff	sequencing	complete_status
AG006	Oryza_sativa	oryza	Oryza_clonal_II	Italy	2018	data/references/AG006.fa	data/annotations/AG006.gff3	ONT_Illumina	complete
BR62	Eleusine_indica	eleusine	Eleusine	NA	NA	data/references/BR62.fa	data/annotations/BR62.gff3	ONT_Illumina	complete
...
```

### 5.2 Test isolate metadata

File: `config/samples.tsv`

```tsv
sample_id	host	host_group	lineage	country	year	platform	fastq	read_n50_bp	estimated_coverage	pilot
TEST01	Oryza_sativa	oryza	Oryza_clonal_I	Italy	2024	ONT	data/longreads/TEST01.fastq.gz	NA	NA	yes
TEST02	Triticum_aestivum	triticum	Triticum	NA	NA	ONT	data/longreads/TEST02.fastq.gz	NA	NA	yes
...
```

### 5.3 Global FASTA headers

All contig names must be unique across all references.

```text
>AG006__chr01
>AG006__chr02
>AG006__mChrA
>BR62__chr01
>BR62__contig07
```

Keep the original FASTA files unchanged. Generate new, renamed working copies and document the mapping in `results/panel/contig_name_map.tsv`.

---

## 6. Phase I: QC and standardization of the 14 references

### 6.1 Assembly QC

```bash
seqkit stats data/references/*.fa > results/qc/assembly_stats.tsv

busco \
  -i data/references/AG006.fa \
  -l sordariomycetes_odb10 \
  -m genome \
  -o AG006_busco \
  -c 16
```

Document per assembly:

- Total size and expected genome size
- Contig and chromosome count
- Contig N50 and largest contigs
- BUSCO completeness, fragmentation, and duplication
- Proportion of unclassified small contigs
- Mitochondrial contigs
- Repeat proportion
- Telomeric repeats at contig ends
- Obvious contamination
- Known or suspected mini-chromosomes

Mitochondria, rDNA arrays, and obvious contamination are not included in the nuclear PAV panel, but are archived separately.

### 6.2 Uniform gene annotation

Use a standardized gene annotation for all 14 genomes wherever possible:

- BRAKER3 or a comparable standardized ab-initio/evidence-based annotation
- Liftoff to transfer high-quality gene models to closely related genomes
- InterProScan or eggNOG-mapper for functional annotation
- SignalP and DeepTMHMM optionally for secreted proteins/effector candidates
- OrthoFinder for orthogroups

Important: liftover alone can underestimate accessory genes. Therefore, supplement it with de novo annotation, especially in non-syntenic, subtelomeric, mini-chromosomal, and repeat-rich regions.

### 6.3 Repeat and TE annotation

Run RepeatModeler2 + RepeatMasker or EDTA consistently on all references. Calculate the repeat proportion per contig and per analysis window.

Repeat-rich regions are not removed, but are flagged as potentially ambiguous. They must not be interpreted as absent solely on the basis of low unique coverage.

---

## 7. Phase II: Pangenomic classification of regions

### 7.1 Orthogroups

```bash
orthofinder \
  -f data/annotations/proteomes/ \
  -S diamond \
  -M msa \
  -T iqtree \
  -t 32 \
  -a 16
```

Use the orthogroup matrix for a gene-based PAV classification. However, a missing gene model is not definitive evidence of true loss. Candidates are therefore additionally validated at the DNA level against the respective assemblies.

### 7.2 Whole-genome alignments

Use multiple anchor references, at least one each associated with Oryza, Triticum, Eleusine, and wild-grass lineages. A single Oryza anchor would systematically represent the diversity of the other lineages as absent or non-syntenic.

```bash
nucmer \
  --maxmatch \
  -l 100 \
  -c 500 \
  -p results/wga/AG006_vs_BR62 \
  data/references/AG006.fa \
  data/references/BR62.fa

delta-filter \
  -m \
  -i 85 \
  -l 500 \
  results/wga/AG006_vs_BR62.delta \
  > results/wga/AG006_vs_BR62.filtered.delta

show-coords \
  -THrd \
  results/wga/AG006_vs_BR62.filtered.delta \
  > results/wga/AG006_vs_BR62.coords.tsv
```

MUMmer4/nucmer is suitable for fast whole-genome alignments; SyRI can distinguish syntenic regions, rearrangements, and local variants on whole-genome alignments. [web:37][web:44]

```bash
syri \
  -c results/wga/AG006_vs_BR62.coords.tsv \
  -r data/references/AG006.fa \
  -q data/references/BR62.fa \
  --prefix results/syri/AG006_vs_BR62_
```

### 7.3 Region types

Each interval in the panel is assigned a class in `results/panel/panel_regions.bed`.

| Region type | Definition |
|---|---|
| `strict_core` | Homologous and largely syntenic in 14/14 reference genomes |
| `soft_core` | Homologous in at least 13/14 reference genomes |
| `shell` | Present in 3–12/14 reference genomes |
| `private_accessory` | Present in only one or two references |
| `accessory_chromosome` | Whole contig/chromosome with low prevalence, lacking stable core synteny, or high structural variability |
| `mini_chromosome` | Small, self-contained chromosome, ideally delimitable by telomeres and clearly separated from the core |
| `subtelomeric_dynamic` | Region near a chromosome end with a high repeat and PAV rate |
| `starship_like` | Mobile candidate block with DUF3435/Captain evidence, cargo gene space, and variable presence |
| `repeat_ambiguous` | Region whose presence is not reliably assessable as binary due to high repeat content or multi-mapping |
| `unclassified` | Not yet sufficiently classified |

Recommended prevalence thresholds for 14 references:

```text
strict core: 14/14
soft core:   13/14
shell:       3–12/14
accessory:   1–2/14
```

### 7.4 Starship-like candidates

A region is only carried as `starship_like` when several lines of evidence converge:

- DUF3435/Captain candidate within or near a putative boundary
- sufficiently large, contiguous region, e.g. at least 20 kb
- cargo genes or functionally heterogeneous gene content
- lacking stable core synteny or strongly variable presence
- optionally terminal repeats or conserved boundaries
- repeat/TE- and/or insertion-rich surroundings

A single DUF3435 hit is not sufficient. Until structural confirmation, the conservative designation is always **Starship-like**.

---

## 8. Phase III: Building the multi-reference panel

### 8.1 Two panel levels

```text
A. Complete catalog:
   14 complete, separate assemblies.

B. Analytical panel:
   - one representative per strict-core homology block
   - alternative soft-core/shell blocks, if not nearly identical
   - all validated accessory regions
   - all mini-chromosomes
   - all Starship-like regions
```

### 8.2 Deduplication rules

- Strict-core homology blocks: choose one representative sequence.
- Soft-core/shell blocks: keep alternative variants if they are biologically or structurally distinct.
- Accessory, mini-chromosome, and Starship-like regions: keep all validated, non-redundant variants.
- Cluster highly identical sequences before inclusion, e.g. at least 98% identity and at least 90% mutual coverage.
- For homologous regions, store all supporting reference isolates as members of a region cluster.

### 8.3 FASTA IDs

```text
>PANEL000104|type=mini_chromosome|cluster=ACC_018|rep=AG006|source=AG006__mChrA
>PANEL000105|type=starship_like|cluster=STAR_004|rep=BR62|source=BR62__chr05:2210000-2460000
>PANEL000106|type=strict_core|cluster=CORE_0321|rep=OryzaRef1|source=OryzaRef1__chr03:100000-450000
```

### 8.4 Panel manifest

File: `results/panel/panel_contig_manifest.tsv`

```tsv
panel_id	region_cluster	region_type	representative_reference	source_interval	length_bp	n_reference_genomes	reference_isolates	host_groups	repeat_fraction	starship_evidence
PANEL000104	ACC_018	mini_chromosome	AG006	AG006__mChrA:1-1200000	1200000	2	AG006;BR62	oryza;eleusine	0.58	no
PANEL000105	STAR_004	starship_like	BR62	BR62__chr05:2210000-2460000	250000	3	BR62;X12;X14	eleusine;wildgrass	0.43	DUF3435+boundary+PAV
PANEL000106	CORE_0321	strict_core	OryzaRef1	OryzaRef1__chr03:100000-450000	350000	14	ALL	all	0.08	no
```

This file is the key to later biological interpretation: it links panel coordinates to region type, original reference coordinate, prevalence, repeat proportion, and all reference isolates carrying a homologous block.

### 8.5 Final panel FASTA

```bash
cat results/panel/panel_core.fa \
    results/panel/panel_accessory.fa \
  > results/panel/Mo_multiref_panel_v1.fa

samtools faidx results/panel/Mo_multiref_panel_v1.fa
```

Version every change to the panel, e.g. `v1.0`, `v1.1`, `v2.0`. Do not retroactively modify old panel versions.

---

## 9. Phase IV: Selection and QC of the ten pilot isolates

### 9.1 Stratified selection

The ten isolates should represent the diversity of the 42 available long-read isolates as well as possible:

- 2–3 Oryza-associated isolates from different clones/populations
- 2 Triticum-associated isolates
- 2 Eleusine-associated isolates
- 2–3 wild-grass or other host isolates
- at least one isolate with expected accessory DNA or mini-chromosomes
- at least one benchmark isolate with an existing complete assembly, if possible

### 9.2 Long-read QC

```bash
NanoPlot \
  --fastq data/longreads/TEST01.fastq.gz \
  --outdir results/qc/TEST01_nanoplot \
  --threads 8

seqkit stats data/longreads/TEST01.fastq.gz \
  > results/qc/TEST01_read_stats.tsv
```

Capture per isolate:

- Total bases and estimated coverage
- Read N50 and mean read length
- Quality distribution
- Proportion of very short reads
- Expected nuclear versus mitochondrial coverage

For robust absence calls, the overall nuclear coverage should be sufficient; a reasonable ONT starting value for the pilot experiment is about 20×, but it must be calibrated based on your data.

---

## 10. Phase V: Mapping the long reads

### 10.1 ONT mapping

```bash
minimap2 \
  -ax map-ont \
  --secondary=yes \
  -t 32 \
  results/panel/Mo_multiref_panel_v1.fa \
  data/longreads/TEST01.fastq.gz \
| samtools sort -@ 8 \
  -o results/mapping/TEST01.panel.bam

samtools index results/mapping/TEST01.panel.bam

samtools flagstat results/mapping/TEST01.panel.bam \
  > results/mapping/TEST01.flagstat.txt

samtools coverage results/mapping/TEST01.panel.bam \
  > results/mapping/TEST01.coverage_by_contig.tsv
```

For PacBio HiFi, `-ax map-hifi` is used.

### 10.2 Two evaluation levels

1. **All alignments including secondary alignments:** evidence of shared homology and evaluation of reference alternatives.
2. **Primary, high-quality, unique alignments:** quantitative PAV evaluation, e.g. MAPQ ≥ 20 or 30.

A region with high homology to multiple references must not be attributed solely on the basis of the primary alignment hit to a single reference.

---

## 11. Phase VI: Window-based PAV analysis

### 11.1 Generating windows

```bash
cut -f1,2 results/panel/Mo_multiref_panel_v1.fa.fai \
  > results/panel/panel.genome

bedtools makewindows \
  -g results/panel/panel.genome \
  -w 10000 \
  > results/panel/panel_10kb_windows.bed
```

Use 10 kb windows as the standard for genome-wide PAV maps. For Starship boundaries, small gene clusters, or breakpoint validation, add a higher resolution with 1–2 kb windows.

### 11.2 Calculating coverage

```bash
mosdepth \
  --by results/panel/panel_10kb_windows.bed \
  --threads 16 \
  results/mapping/TEST01 \
  results/mapping/TEST01.panel.bam
```

Per window, the following are calculated:

- mean and median depth
- breadth at at least 1×, 3×, and 5×
- proportion of high-quality alignments
- proportion of primary versus secondary alignments
- repeat proportion
- number of unique windows or diagnostic k-mers
- number of reads spanning window boundaries

### 11.3 PAV rules

| Region type | Present | Absent | Additional rule |
|---|---|---|---|
| Strict-/soft-core | Breadth ≥ 0.80 at at least 5× | Breadth < 0.10 | Only with sufficient global coverage |
| Accessory block | Breadth ≥ 0.70 in at least 80% of unique windows | Breadth < 0.10 in at least 80% of unique windows | Evaluate multi-mapping explicitly |
| Mini-chromosome | Support in multiple unique windows distributed along the length | No stable coverage across length and clear absence of markers | Do not derive from individual TE windows |
| Starship-like | Interior region plus at least one boundary supported | Interior region and boundaries not supported | Boundary support increases confidence |
| Repeat-ambiguous | No binary call | No binary call | Status `ambiguous` |

Adjacent windows with the same status are merged into PAV blocks. A project-wide minimum block length can initially be set at 20 kb and later adjusted based on benchmark data.

### 11.4 Evidence classes

```text
high_confidence:
  Coverage across unique windows + high-quality alignments +
  breakpoint evidence or multiple spanning reads

moderate_confidence:
  Stable coverage across unique windows, but without breakpoint evidence

ambiguous:
  strong repeat/multi-mapping signals, conflicting evidence, or
  insufficient global long-read depth
```

---

## 12. Phase VII: SV calling with long reads

In addition to coverage, breakpoint-based SV calling is performed per test isolate.

```bash
sniffles \
  --input results/mapping/TEST01.panel.bam \
  --vcf results/pav/TEST01.sv.vcf.gz \
  --snf results/pav/TEST01.snf \
  --reference results/panel/Mo_multiref_panel_v1.fa \
  --threads 24
```

For the cohort:

```bash
sniffles \
  --input results/pav/TEST01.snf \
          results/pav/TEST02.snf \
          results/pav/TEST03.snf \
          results/pav/TEST04.snf \
          results/pav/TEST05.snf \
          results/pav/TEST06.snf \
          results/pav/TEST07.snf \
          results/pav/TEST08.snf \
          results/pav/TEST09.snf \
          results/pav/TEST10.snf \
  --vcf results/pav/pilot10.cohort.sv.vcf.gz \
  --threads 24
```

Sniffles2 can detect long-read-based deletions, insertions, duplications, inversions, and translocations, and allows population-wide merging via `.snf` files. [web:20][web:22]

The final PAV decision combines:

```text
PAV evidence = coverage breadth + mapping uniqueness + SV breakpoints + spanning reads
```

---

## 13. Phase VIII: Registering and assigning candidate regions

### 13.1 Defining candidate regions before mapping

Each candidate region receives a stable ID:

```text
ACC_0001 ... ACC_n       accessory region
MCHR_0001 ... MCHR_n     mini-chromosome or mChr block
STAR_0001 ... STAR_n     Starship-like region
SUBTEL_0001 ... SUBTEL_n subtelomeric dynamic region
```

File: `results/panel/panel_candidate_regions.tsv`

```tsv
candidate_id	panel_id	region_type	panel_start	panel_end	length_bp	representative_reference	reference_isolates	n_reference_genomes	host_groups	core_accessory_class	repeat_fraction	captain_gene_id	cargo_gene_count	annotation_confidence
MCHR_0001	PANEL000104	mini_chromosome	1	1200000	1200000	AG006	AG006;BR62	2	oryza;eleusine	accessory_chromosome	0.58	NA	42	high
STAR_0004	PANEL000105	starship_like	1	250000	250000	BR62	BR62;X12;X14	3	eleusine;wildgrass	starship_like	0.43	BR62_g08765	17	medium
ACC_0017	PANEL000221	accessory_region	1	85000	85000	OryzaRef3	OryzaRef3;AG006	2	oryza	private_accessory	0.31	NA	6	high
```

### 13.2 Result per test isolate

File: `results/pav/candidate_region_calls.tsv`

```tsv
test_isolate	candidate_id	region_type	status	confidence	breadth_5x	mean_depth	unique_window_fraction	sv_support	shared_reference_isolates	closest_panel_reference	assignment_confidence	notes
TEST01	MCHR_0001	mini_chromosome	present	high	0.94	31.2	0.89	12	AG006;BR62	AG006	medium	mChrA-like region
TEST01	STAR_0004	starship_like	absent	high	0.03	0.4	0.96	0	BR62;X12;X14	NA	NA	no interior or boundary support
TEST02	ACC_0017	accessory_region	ambiguous	low	0.43	6.1	0.28	1	OryzaRef3;AG006	OryzaRef3	low	repeat-associated multi-mapping
```

The columns `shared_reference_isolates` and `closest_panel_reference` must be distinguished:

- `shared_reference_isolates`: all reference isolates carrying the homologous region cluster.
- `closest_panel_reference`: the reference variant with the best alignment evidence.
- `assignment_confidence`: confidence of this specific assignment.

For nearly identical regions, "mChrA-like and shared in AG006/BR62" is more robust than claiming a definitive, donor-specific origin.

### 13.3 Localization and biological annotation

```bash
bedtools intersect \
  -a results/pav/TEST01.pav_blocks.bed \
  -b results/panel/panel_regions.bed \
  -wa -wb \
  > results/pav/TEST01.pav_blocks_localized.tsv

bedtools intersect \
  -a results/pav/TEST01.pav_blocks.bed \
  -b results/panel/panel_genes.gff3 \
  -wa -wb \
  > results/pav/TEST01.pav_gene_overlap.tsv
```

Every final PAV call should contain:

```text
PAV_ID
Test isolate
Panel coordinate
Length
PAV status
Confidence class
Region type
Core/accessory class
Subtelomeric yes/no
Repeat proportion
Overlapping genes
Orthogroups
Functional annotation
Starship/Captain evidence
Reference cluster
Shared reference isolates
Closest panel reference
```

---

## 14. Phase IX: Rarefaction and panel saturation

Rarefaction assesses two distinct questions:

1. How many new candidate regions are discovered with additional reference genomes?
2. Does the reference panel adequately cover the candidate regions of the long-read test isolates?

### 14.1 Reference rarefaction

For each panel size \(k = 1, ..., 14\):

1. Draw many random combinations of \(k\) reference genomes.
2. Merge the candidate region clusters present in the combination.
3. Count non-redundant regions separately for accessory, mini-chromosome, and Starship-like regions.
4. Calculate mean, median, standard deviation, and 95% interval.

Formally:

\[
R(k) = \left| \bigcup_{i \in S_k} C_i \right|
\]

Here, \(C_i\) is the set of candidate regions of reference \(i\), \(S_k\) is a selection of \(k\) references, and \(R(k)\) is the number of accumulated, non-redundant candidate regions discovered.

Proposal for `config/thresholds.yaml`:

```yaml
rarefaction:
  n_permutations: 1000
  saturation_delta_fraction: 0.05
```

Expected output: `results/rarefaction/rarefaction_reference_panel.tsv`

```tsv
panel_size	region_class	mean_n_regions	median_n_regions	ci_lower	ci_upper	n_permutations
1	accessory	24	24	18	30	1000
1	starship_like	3	3	1	5	1000
...
14	accessory	218	218	218	218	1
14	starship_like	31	31	31	31	1
```

### 14.2 Saturation criterion

Calculate the additional gain of the last panel size:

\[
\Delta R_{13 \rightarrow 14} =
\frac{R(14)-R(13)}{R(14)}
\]

A panel is assessed as approximately saturated if:

- the mean gain from adding the last references is below approximately 2–5%,
- the confidence intervals become narrow,
- and the long-read test isolates only rarely yield new, longer, and well-supported candidate sequences without panel homology.

Report separate curves for:

- accessory regions overall
- mini-chromosomes
- Starship-like regions
- subtelomeric dynamic regions
- all candidate regions combined

A plateau in the overall curve does not automatically mean saturation for mini-chromosomes or Starship-like regions; these classes in particular can accumulate considerably more slowly.

### 14.3 Test-isolate rarefaction and novelty check

For test isolates with sufficient sequencing depth:

1. Extract unmapped or poorly mapped reads.
2. Perform targeted local assemblies.
3. Search the resulting contigs against the multi-reference panel.
4. Define well-supported, panel-external sequences as new candidate regions.

A novel candidate region could, for example, fulfill:

- local contig at least 10–20 kb
- not predominantly repetitive sequence
- no panel homology over at least 80% of the length
- accessory, mini-chromosome, Starship-like, or subtelomeric candidate property

Then generate a second accumulation curve:

```text
x-axis: number of test isolates analyzed
y-axis: cumulative number of new candidate regions not represented in the panel
```

If this curve continues to rise sharply, the panel is incomplete. If it flattens early, this indicates good coverage of the relevant candidate space.

---

## 15. Configurable parameters

File: `config/thresholds.yaml`

```yaml
panel:
  strict_core_fraction: 1.00
  soft_core_fraction: 0.93
  shell_min_fraction: 0.21
  candidate_min_length_bp: 10000
  dedup_identity: 0.98
  dedup_coverage: 0.90

mapping:
  preset_ont: map-ont
  preset_hifi: map-hifi
  min_mapq_unique: 20
  min_primary_alignment_bp: 3000

pav:
  window_size_bp: 10000
  high_resolution_window_size_bp: 2000
  min_depth_present: 5
  present_breadth_core: 0.80
  present_breadth_accessory: 0.70
  absent_breadth: 0.10
  min_block_bp: 20000
  min_global_depth: 20

starship:
  min_region_length_bp: 20000
  duf3435_evalue: 1.0e-5
  require_cargo_genes: true
  require_boundary_evidence_for_high_confidence: true

rarefaction:
  n_permutations: 1000
  saturation_delta_fraction: 0.05
```

All thresholds are starting values and must be calibrated with benchmark data. Changes to parameters create a new analysis version and are documented in the provenance file.

---

## 16. Validation before scaling to 42 isolates

Before the full analysis of the remaining 32 long-read isolates:

1. Use at least one to several isolates with a known high-quality assembly as a benchmark.
2. Temporarily remove the respective assembly from the reference panel.
3. Map the long reads of this isolate against the remaining panel.
4. Compare the PAV calls with assembly-vs-assembly analyses.
5. Calculate sensitivity, precision, and F1 separately for core, accessory, mini-chromosomal, subtelomeric, and repeat-rich regions.
6. Check the stability of the findings under moderate parameter variation, e.g. MAPQ 20 versus 30 and breadth 0.70 versus 0.80.

Initial target values can be:

- sensitivity of at least 90% for clearly defined PAVs above a minimum size
- specificity of at least 95% in unambiguous, non-repetitive regions
- conservative handling of repeats and subtelomeric regions

---

## 17. Cost and compute strategy

| Phase | Effort | Approach |
|---|---:|---|
| FASTA QC, headers, metadata | low | Carry out immediately for all 14 references |
| Annotation and repeat masking | medium to high | Once, consistently, for all references |
| Whole-genome alignments/SyRI | medium | Once per sensible reference combination |
| Panel and candidate construction | medium | Versioned; repeat only on panel update |
| Mapping of the 10 pilot isolates | medium to high | Pilot before the complete 42-isolate cohort |
| Coverage-based PAV | low to medium | For all pilot isolates |
| SV calling | medium | After successful mapping and QC |
| Local assemblies of unmapped reads | high | Only targeted, for novelty/HGT/PAV candidates |
| Scaling to 42 isolates | scalable | Only after validation of the pilot parameters |

---

## 18. Decision logic for final findings

A biological statement about a candidate region should fulfill the following chain:

```text
Panel region defined
→ Region class and reference cluster known
→ Sufficient long-read quality in the test isolate
→ Coverage across unique windows
→ Assessment of multi-mapping and repeats
→ Optional SV/breakpoint or spanning-read evidence
→ PAV status and confidence class
→ Assignment to all sharing reference isolates
```

### Interpretation

| Finding | Permissible interpretation |
|---|---|
| High-confidence presence of a mini-chromosome cluster | The test isolate carries an mChr-like region homologous to the reference isolates listed in the cluster |
| Presence of a Starship-like cluster including boundary support | The test isolate likely carries a homologous Starship-like region; structural confirmation may require additional local assembly |
| Absence of a unique accessory region with sufficient depth | Robust loss or non-presence in the test isolate |
| Low coverage in TE-rich, multiply homologous DNA | Ambiguity; no firm absence statement |
| Best hit on one reference isolate, but multiple references in the cluster | Closest reference variant, but no definitive donor or origin assignment |

---

## 19. Concluding principle

The multi-reference panel is not just a FASTA file. It is a versioned catalog of homologous sequence regions with:

- stable region and cluster IDs,
- coordinates in the analytical panel,
- original coordinates in the 14 reference genomes,
- core/accessory/mini-chromosome/Starship-like classification,
- repeat and gene annotation,
- prevalence in the reference set,
- and a complete list of all reference isolates sharing the region.

This allows every long-read PAV to later be traceably localized as core, accessory, mini-chromosomal, subtelomeric, or Starship-like. At the same time, it allows transparent documentation of which reference isolates also carry a candidate region found in a test isolate, and whether the panel is already sufficiently saturated for the diversity of the blast lineages studied.
