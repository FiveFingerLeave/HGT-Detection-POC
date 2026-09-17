# Development Environment and AI Integration for the Barragan Pipeline

**Project:** PhD thesis – Barragan / *Magnaporthe oryzae*  
**Date:** August 31, 2026  
**Purpose:** Building a reproducible working environment for Starfish, Starship, mini-chromosome, and isolate comparison analyses.

---

## 1. Goal of the Setup

For this PhD project, the analysis should not consist of individually copied terminal commands. The goal is a version-controlled, reproducible bioinformatics project that can be scaled from a single test isolate to many isolates.

The core technical idea is:

```text
VS Code as the working interface
        +
Git for version control
        +
Conda/Mamba for software environments
        +
Snakemake for automated pipelines
        +
Python/R for custom data processing and analysis
        +
AI as controlled support for code, tests, and documentation
```

AI can help with explanations, code drafts, debugging, and documentation. However, the actual scientific responsibility remains with the project: input data, versions, parameters, tests, and biological interpretation must be traceably recorded and verified.

---

## 2. Recommended Overall Architecture

```text
Work computer / browser
        │
        ├── Visual Studio Code
        │   ├── Editor for Python, Snakemake, YAML, and Markdown
        │   ├── Integrated Linux terminal
        │   ├── Git view and commit history
        │   └── AI chat / pair programming
        │
        └── Remote-SSH, if the compute environment is a server
                    │
                    ▼
Linux compute environment: MilleniumFalke
        │
        ├── Project repository: ~/promotion/barragan/
        ├── Conda/Mamba environments
        ├── Snakemake workflow
        ├── Raw data outside Git
        ├── Result data outside Git
        └── Logs, benchmarks, and reports
```

If `MilleniumFalke` is a separate Linux server, VS Code is installed locally on the work computer and connected via SSH. If VS Code runs directly on `MilleniumFalke`, the project folder can be opened directly and locally.

---

## 3. Visual Studio Code

### Why VS Code?

VS Code is a suitable working interface because it combines editor, terminal, Git support, debugging, and extensions for Python, YAML, Markdown, and remote development.

The following file types are primarily edited in this project:

| File type | Use |
|---|---|
| `.py` | small scripts for formatting, validation, tables, and analysis |
| `Snakefile`, `.smk` | Snakemake pipeline rules |
| `.yaml`, `.yml` | parameters, software environments, and reference paths |
| `.tsv` | sample metadata, mapping tables, result summaries |
| `.md` | documentation, decisions, method drafts, and lab notebook |
| `.sh` | individual helper commands, where more sensible than Python |

### Recommended Extensions

| Extension | Purpose |
|---|---|
| Remote - SSH | Opens a folder on a Linux machine/server via SSH |
| Python | Python code, interpreter selection, linting, and debugging |
| Pylance | Autocompletion and type checking for Python |
| Snakemake | Syntax highlighting and support for Snakemake files |
| YAML | Syntax checking and editing of YAML configurations |
| Markdown All in One | Markdown preview, outline, and formatting |
| GitLens | Extended Git history and file versions |
| Error Lens | Shows errors and warnings directly in the editor |
| Jupyter | Optional for exploratory notebooks; not for final pipeline steps |

### Opening the Project Folder

Directly on the Linux system:

```bash
code ~/promotion/barragan
```

If the `code` command is not yet available, the folder

```text
~/promotion/barragan
```

can be opened in VS Code via `File → Open Folder…`.

---

## 4. Setting Up Remote-SSH

If `MilleniumFalke` is reached via SSH, Remote-SSH is the most convenient way to work. The editor then operates on the server: files, terminal commands, Conda environments, and extensions use the Linux compute environment directly.

### Creating an SSH Alias

On the local work computer, the file `~/.ssh/config` can contain an alias:

```sshconfig
Host milleniumfalke
    HostName YOUR_SERVER_OR_YOUR_IP
    User flori
    IdentityFile ~/.ssh/id_ed25519
    ServerAliveInterval 60
```

The connection then works in the terminal with:

```bash
ssh milleniumfalke
```

In VS Code:

1. Open the Command Palette with `Ctrl+Shift+P`.
2. Select `Remote-SSH: Connect to Host...`.
3. Select the host `milleniumfalke`.
4. After a successful connection, open `~/promotion/barragan` as the folder.

---

## 5. Conda and Mamba

### Role of Conda/Mamba

Bioinformatics tools often have many dependencies. Conda and Mamba ensure that defined tool versions can be installed together.

- **Conda** manages software environments.
- **Mamba** is largely compatible with Conda but often resolves dependencies faster.
- **One environment per task class** prevents conflicts between programs.

### Existing and Planned Environments

| Environment | Task |
|---|---|
| `starfish_env` | existing environment for Starfish, MetaEuk, HMMER, and associated databases |
| `barragan-workflow` | Snakemake, Python, and pipeline control |
| `alignment` | later e.g. minimap2, samtools, MUMmer/nucmer, Mash/Skani |
| `barragan-python` | table analysis, validation, custom scripts, and figures |

### Checking Availability

```bash
git --version
conda --version
mamba --version
snakemake --version
python3 --version
```

If Mamba is missing, it can be installed in the base environment:

```bash
conda install -n base -c conda-forge mamba
```

### Creating the Workflow Environment

```bash
mamba create -n barragan-workflow \
  -c conda-forge -c bioconda \
  snakemake \
  python=3.12 \
  pandas \
  biopython \
  pyyaml \
  pytest \
  ruff \
  git \
  -y
```

Activating:

```bash
conda activate barragan-workflow
```

Starfish initially remains in `starfish_env`. In the long term, the pipeline defines its own environment files so that Snakemake can provide the correct environment per rule.

---

## 6. Git Repository

### Why Git?

Git stores the history of the code and configurations. This makes it possible to answer:

- When was a parameter changed?
- Which version of a script produced a particular result?
- Why was a rule added or removed?
- How can a working state be restored?

Git stores **code, configuration, and documentation**, but generally not large FASTA, FASTQ, BAM, or complete result files.

### Initializing the Repository

```bash
cd ~/promotion/barragan
git init
```

### Recommended Structure

```text
barragan/
├── README.md
├── .gitignore
├── .github/
│   └── copilot-instructions.md
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
│   │   ├── minichromosomes.smk
│   │   ├── comparison.smk
│   │   └── report.smk
│   ├── scripts/
│   │   ├── normalize_fasta_headers.py
│   │   ├── normalize_gff_seqids.py
│   │   ├── validate_fasta_gff_ids.py
│   │   ├── calculate_contig_metrics.py
│   │   ├── parse_starfish_results.py
│   │   ├── build_pav_matrix.py
│   │   └── score_hgt_candidates.py
│   └── envs/
│       ├── starfish.yaml
│       ├── alignment.yaml
│       ├── python.yaml
│       └── qc.yaml
├── input/
│   ├── manifests/
│   └── metadata/
├── results/
│   ├── qc/
│   ├── normalized/
│   ├── starfish/
│   ├── minichromosomes/
│   ├── comparisons/
│   └── reports/
├── logs/
├── benchmarks/
├── docs/
│   ├── workflow.md
│   ├── decisions.md
│   ├── data_dictionary.md
│   ├── methods_draft.md
│   └── ai_usage.md
└── tests/
    ├── data/
    └── test_normalization.py
```

### `.gitignore`

Create a `.gitignore` as early as possible:

```gitignore
# Large raw data and sequence data
assemblies/
raw_data/
reads/
*.fastq
*.fastq.gz
*.fq
*.fq.gz
*.bam
*.bai
*.cram
*.crai

# Large bioinformatics files
*.fna
*.fa
*.fasta
*.fas
*.faa
*.gff
*.gff3
*.hmm
*.dmnd
*.mmi

# Automatically generated result and temporary files
results/
logs/
benchmarks/
temp/
tmp/
.snakemake/

# Python and editor artifacts
__pycache__/
*.pyc
.ipynb_checkpoints/
.vscode/settings.json
.DS_Store

# Credentials and secrets
.env
*.key
```

### First Commit

```bash
git add README.md .gitignore config workflow docs tests .github
git commit -m "Initialize reproducible Starfish and minichromosome workflow"
```

Afterward, commit in small, topically comprehensible units:

```bash
git status
git add workflow/scripts/normalize_fasta_headers.py tests/test_normalization.py
git commit -m "Add FASTA header normalization with ID mapping"
```

---

## 7. Snakemake as the Pipeline Engine

### Role of Snakemake

Snakemake describes which result files are produced from which inputs. It recognizes dependencies between steps, can compute multiple isolates in parallel, and re-runs missing or incomplete steps.

A simplified pipeline:

```text
Original FASTA
      ↓
Normalize FASTA headers
      ↓
Synchronize GFF SeqIDs
      ↓
Validate FASTA/GFF compatibility
      ↓
Starfish: search for YR candidates
      ↓
Analyze Starship context
      ↓
Determine contig metrics and mChr candidates
      ↓
Isolate comparison: PAV, homology, core-vs-element contrast
```

### Starting Small

The first pilot workflow should contain only four steps:

1. Normalize FASTA headers.
2. Synchronize GFF SeqIDs with the same ID mapping table.
3. Check whether every GFF SeqID occurs in the matching FASTA.
4. Start the Starfish YR annotation.

Only once this process works stably for one isolate do three isolates follow, and then the full dataset.

### Example: Central Snakefile

File `workflow/Snakefile`:

```python
configfile: "../config/parameters.yaml"

include: "rules/headers.smk"
include: "rules/starfish.smk"

rule all:
    input:
        expand(
            "../results/starfish/{sample}/{sample}_YR.filt.gff",
            sample=config["pilot_samples"],
        )
```

### Example: Parameter File

File `config/parameters.yaml`:

```yaml
project_name: barragan_hgt
separator: "_"

threads:
  starfish: 6

pilot_samples:
  - GCA_004346965.1

starfish:
  profile: /home/flori/miniforge3/envs/starfish_env/db/YRsuperfams.p1-512.hmm
  proteins: /home/flori/miniforge3/envs/starfish_env/db/YRsuperfamRefs.faa
  hmm_evalue: 0.001
```

Absolute paths are acceptable for the first pilot run but must be documented. Later, reference files should be managed in a controlled resource structure or via central path configurations.

### Dry Run Before Real Runs

Before a computationally intensive workflow is started:

```bash
cd ~/promotion/barragan
conda activate barragan-workflow

snakemake \
  --snakefile workflow/Snakefile \
  --cores 6 \
  --use-conda \
  --dry-run \
  --printshellcmds
```

A dry run shows planned commands and missing inputs without generating files.

A real run:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --cores 6 \
  --use-conda \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds
```

| Option | Meaning |
|---|---|
| `--cores 6` | use a maximum of six compute cores |
| `--use-conda` | use the declared software environment per rule |
| `--rerun-incomplete` | cleanly regenerate incomplete results |
| `--keep-going` | continue computing independent isolates if one fails |
| `--printshellcmds` | make executed terminal commands visible |

---

## 8. AI Integration

### Suitable AI Tasks

AI is particularly helpful for:

- Explaining error messages, tool help, and file formats
- Drafting small Python scripts
- Converting a tested terminal command into a Snakemake rule
- Writing and improving unit tests
- Code review and searching for edge cases
- Creating YAML, TSV, and Markdown templates
- Formulating method documentation
- Detecting missing input/output/log declarations in workflow rules

### Tasks That Must Always Be Checked Personally

AI must not decide unchecked on:

- biological interpretation of individual hits
- final HGT or Starship calls
- thresholds without calibration
- deleting, renaming, or overwriting large data holdings
- tool parameters without local help, literature, or a test run
- private raw data, credentials, API keys, or unpublished sensitive metadata

### Safe Order of Work

```text
AI suggestion
      ↓
Check local tool help / primary source
      ↓
Test with a small dataset
      ↓
Validate output, logs, and edge cases
      ↓
Test and review code
      ↓
Commit the change
      ↓
Only then scale to many isolates
```

### GitHub Copilot in VS Code

When using GitHub Copilot, a project-specific rules file can be created:

```bash
mkdir -p .github
code .github/copilot-instructions.md
```

Recommended content:

```markdown
# Barragan bioinformatics workflow instructions

## Scientific and reproducibility rules
- Never alter raw input files in place.
- Write derived files only to results/ or temporary directories.
- Every rule must declare explicit inputs, outputs, logs, threads, and conda environment.
- Do not infer biological conclusions from a single YR hit.
- Treat YR hits as YR_candidate unless additional evidence is documented.
- Preserve isolate IDs and record every ID transformation in a mapping table.
- For FASTA header changes, generate old-to-new ID mappings.
- For GFF changes, validate that every sequence ID exists in the corresponding FASTA.
- Do not hard-code sample names in scripts; read them from config or samples.tsv.
- Do not use destructive commands without first proposing a safe alternative.
- Add or update tests for each custom Python transformation.
- Use clear docstrings, type hints, TSV outputs, and deterministic sorting.

## Code style
- Use Python 3.12 syntax.
- Use pathlib instead of shell-specific path construction where possible.
- Use pandas only for tabular data; use Biopython for FASTA parsing.
- Fail clearly with actionable error messages.
- Keep scripts small: one biological/data transformation per script.
```

This file is a guideline for the AI. It is meant to prevent it from altering raw data, hard-coding samples into the code, or labeling preliminary YR hits as confirmed Starships.

### Example Prompts

**Python script for FASTA headers:**

```text
Write a Python 3.12 script for workflow/scripts/normalize_fasta_headers.py.
Input: FASTA and isolate_id.
Output: normalized FASTA with headers <isolate_id>_<original_id> and a TSV mapping file old_id, new_id.
Do not alter sequences. Abort with an understandable error message if a header is empty or duplicated.
Also create pytest tests with a minimal test FASTA.
```

**Snakemake rule:**

```text
Create a Snakemake rule normalize_fasta_headers.
The input comes from config/samples.tsv, output is results/normalized/{sample}.fna and results/normalized/{sample}.id_map.tsv.
Write a log to logs/headers/{sample}.log.
Use workflow/envs/python.yaml, declared threads, and no hard-coded isolate IDs.
Then briefly explain each line.
```

**Review of a data format:**

```text
Critically review this GFF and FASTA ID mapping logic.
List possible failure cases: FASTA headers with spaces, duplicate contig IDs, GFF SeqIDs without a FASTA counterpart, comments in the GFF, and compressed input files.
Propose concrete tests, but do not change any files.
```

---

## 9. Documentation and Digital Lab Notebook

In addition to the pipeline, the project needs a traceable methods and decisions log.

### `README.md`

The README contains:

- Goal of the pipeline
- Prerequisites and installation
- Brief folder overview
- Example of a test run
- Most important outputs

### `docs/workflow.md`

Documents:

- Order of all pipeline rules
- Input and output of each rule
- Tools used and their purpose
- Expected resources and runtimes

### `docs/decisions.md`

Documents methodological decisions, e.g.:

```markdown
## 2026-08-31 — Unique Contig IDs

**Decision:** FASTA headers are normalized to
`<assembly_accession>_<original_contig_id>`.

**Rationale:** Starfish expects a genome ID plus a feature/contig ID.
Unique IDs prevent collisions between isolates.

**Consequence:** Corresponding GFF SeqIDs are synchronously adjusted via the
same mapping table and subsequently validated against the FASTA.
```

Additionally document:

- Starfish version
- Starfish database files
- HMM thresholds
- MetaEuk parameters
- Criteria for mChr classes
- Criteria for Starship confidence
- Criteria for shared elements and HGT prioritization

### `docs/ai_usage.md`

AI usage is transparently recorded here:

```markdown
# AI Usage Log

## Principle
AI serves for explanation, code drafts, documentation, and review.
All scientific decisions, parameters, tests, and interpretations
are checked and documented by me.

## 2026-08-31
- Use: Explanation of the Starfish workflow and drafting of a
  documentation structure.
- Verification: Local `starfish annotate --help` output checked;
  single-isolate test run executed.
- Result: 13 HMM-validated YR candidates in GCA_004346965.1.
- Interpretation: Documented as YR_candidate, not as a confirmed Starship.
```

### Data Provenance and Checksums

For each assembly, accession, origin, download date, file, GFF version, and checksum should be stored.

Example:

```bash
sha256sum \
  assemblies/ncbi_dataset/ncbi_dataset/data/GCA_004346965.1/GCA_004346965.1_ASM434696v1_genomic.fna \
  > input/manifests/GCA_004346965.1.sha256
```

---

## 10. Concrete Starting Plan

### Today: Create the Basic Structure

```bash
cd ~/promotion/barragan

mkdir -p \
  config \
  workflow/rules \
  workflow/scripts \
  workflow/envs \
  input/manifests \
  input/metadata \
  results \
  logs \
  benchmarks \
  docs \
  tests/data \
  .github

touch \
  README.md \
  .gitignore \
  config/samples.tsv \
  config/parameters.yaml \
  workflow/Snakefile \
  docs/decisions.md \
  docs/ai_usage.md \
  .github/copilot-instructions.md

git init
```

Then open the project folder in VS Code:

```bash
code ~/promotion/barragan
```

### This Week: A Reproducible Single-Isolate Pilot

1. Complete `config/samples.tsv` for `GCA_004346965.1`.
2. Write and test the Python script for FASTA header normalization.
3. Generate the mapping table `old_id → new_id`.
4. Find the matching NCBI GFF and adjust its SeqIDs with the same mapping table.
5. Implement FASTA/GFF validation.
6. Repeat Starfish with the normalized inputs.
7. Document log, parameters, result, and interpretation as `YR_candidate`.
8. Create a Git commit after each working stage.

### Afterward: Three-Isolate Pilot

1. Add two contrasting additional isolates to `samples.tsv`.
2. Run the workflow without manual path changes.
3. Merge YR candidates into an overall table.
4. Calculate contig lengths, GC content, and gene counts as initial mChr metrics.
5. Test initial pairwise alignments and PAV logic with the three isolates.
6. Scale to the full dataset only after technical and biological quality control.

---

## 11. What Is Not Needed for Now

Not required to get started:

- Docker or Kubernetes
- A complex database solution
- A comprehensive machine learning stack
- AI agents with write access to all data folders
- Fully automatic HGT decisions

For now, clear TSV/YAML files, a small Snakemake project, tested Python scripts, and consistent documentation are sufficient.

---

## 12. Most Important First Milestone

The first technical success is not analyzing all isolates simultaneously. The success is:

> One isolate runs automatically from the original FASTA through normalized, validated IDs to a traceable Starfish YR output via Snakemake.

Only once this process is stable, tested, and documented should the pipeline be extended to multiple isolates and then to the complete dataset.
