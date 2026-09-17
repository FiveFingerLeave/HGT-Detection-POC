# Documentation: Preparation of the POC for HGT and Starship Analysis in *Magnaporthe oryzae*

## Documentation Status

**Project:** PhD thesis – Population genomics of horizontal gene transfer (HGT) in *Magnaporthe oryzae*  
**Work phase:** Building a global, assembly-based candidate pool for a proof of concept (POC)  
**Working environment:** Windows Terminal with Ubuntu via WSL; Miniforge/conda; isolated conda environment `starfish_env`  
**Documentation status:** Prior to downloading and structurally analyzing the selected genomes

---

## 1. Objective

This PhD project investigates horizontal gene transfer (HGT) in *Magnaporthe oryzae* and *M. grisea*. It focuses on two classes of large mobile genomic elements:

- accessory mini-chromosomes (mChr)
- Starships and starship-like large mobile elements

The overarching hypothesis is that HGT does not occur purely at random, but arises from recurring relationships between donor lineages, recipient lineages, mobile vectors, and local ecological conditions. The subsequent analysis aims to identify global donor–recipient–vector relationships, test their local recurrence in Italian field populations, and model HGT hotspots using generalized linear mixed models (GLMM).

The current work step is explicitly **not a complete HGT demonstration and not a complete global analysis**. It is a preparatory proof of concept that first addresses the following feasibility question:

> Are there enough publicly available, sufficiently contiguous, and ideally long-read-based *M. oryzae* assemblies to reliably detect Starships and mChr structures and to build a representative structural reference panel from them?

Answering this question is important because Starship boundaries, Captain genes, flanking sequences, insertion sites, and mChr contexts often cannot be reliably reconstructed in heavily fragmented short-read assemblies.

---

## 2. Methodological Framework

### 2.1 Separation of Two Levels of Analysis

The workflow deliberately separates a structural assembly level from a global prevalence level.

| Level of analysis | Purpose | Primary data basis | Suitable methods |
|---|---|---|---|
| Structural discovery | Identify complete or largely complete Starships, mChr candidates, Captains, boundaries, flanking regions, and structural variants | Long-read or hybrid assemblies with high contiguity | Starfish/Stargraph, sequence comparison, synteny/alignments, mChr characterization |
| Global prevalence | Measure the occurrence, presence/absence, and cross-lineage distribution of the previously defined candidates across the broad isolate panel | Short reads or existing mapping/coverage data | Mapping against references, breadth-of-coverage, PAV analysis, SNP/haplotype comparison |

This separation prevents a methodologically problematic conflation of "partially visible in a fragmented assembly" with "structurally documented as a complete mobile element."

### 2.2 POC Logic

The POC follows a decision logic with four modules:

1. **Global frequency POC:** Identification of possible identical or highly similar mChr/Starship elements across different lineages.
2. **Local recovery approximation:** Estimation of how many of the global candidates are expected to be recovered in local Italian field populations.
3. **GLMM power analysis:** Simulation of the number of true HGT events required for a stable hotspot model.
4. **Alternative analysis strategies:** Simplifying the model, pooling vector classes, or exploratory/descriptive analysis if the number of events proves insufficient.

The current status concerns exclusively the technical and data-related preparation of Module 1.

---

## 3. Data Basis

### 3.1 Barragán Supplementary Data

The project-specific starting file is `Copy of Barragan2024_SupplementalTables.xlsx`. The relevant sheets are:

| Sheet | Content | Use in the current workflow |
|---|---|---|
| `TableS1` | Metadata of nine Italian clonal rice blast isolates | Local context: isolate, host, location, year, and coordinates |
| `TableS2` | Metadata of 274 isolates with genome-wide SNPs | Global background and strain/host metadata |
| `TableS3` | Contig counts and contig lengths of nine Italian assemblies | Initial assessment of local assembly fragmentation and mChr contigs |
| `TableS4` | Summary of the mChr contigs in the nine Italian isolates | mChr context and known candidates |
| `TableS5` | Metadata of 413 *M. oryzae* and *M. grisea* isolates | Central reference table for later linkage: strain ID, host, lineage, origin, BioProject, and BioSample |
| `TableS6` | Genome-wide and mChr-specific breadth-of-coverage values across 413 isolates | Basis for the later global coverage-/PAV-based analysis |
| `TableS10` | Long-read sequencing and assembly statistics for the *Eleusine* isolate Br62 | Documented long-read reference with quality metrics |

`TableS5` in particular contains the columns `Strain ID`, `Host plant`, `M. oryzae Lineage`, `BioProject`, and `Sample`. The `Sample` column contains BioSample IDs such as `SAMN...` or `SAMEA...`. These are **not assembly accessions**, but identifiers for biological samples and metadata.

### 3.2 NCBI Assembly Data

No FASTA files were downloaded for the POC. Instead, only NCBI assembly metadata for the taxon `Magnaporthe oryzae` were queried initially.

The following was executed:

```bash
./datasets summary genome taxon "Magnaporthe oryzae" \
  --as-json-lines \
  > ncbi_m_oryzae_assemblies.jsonl
```

Result:

- **605 public NCBI assembly records**
- File: `ncbi_m_oryzae_assemblies.jsonl`
- File size: approximately **1.5 MB**
- Only metadata were loaded; no genomes, annotations, or raw reads.

The subsequent TSV conversion captured the following fields:

```bash
./dataformat tsv genome \
  --inputfile ncbi_m_oryzae_assemblies.jsonl \
  --fields accession,organism-name,assminfo-name,assminfo-level,assminfo-assembly-method,assminfo-sequencing-tech,assminfo-bioproject,assminfo-biosample-accession,assminfo-release-date,assminfo-submitter,assmstats-total-sequence-len,assmstats-number-of-contigs,assmstats-contig-n50,assmstats-scaffold-n50,assmstats-genome-coverage \
  > ncbi_m_oryzae_assemblies.tsv
```

This produced the file `ncbi_m_oryzae_assemblies.tsv`. It contains one row per NCBI assembly and forms the basis for the automated selection.

---

## 4. Identifiers and Data Formats Used

| Prefix/format | Meaning | Use |
|---|---|---|
| `GCA_...` | GenBank assembly accession | Unique identifier for a concretely available genome assembly; correct input for assembly download |
| `GCF_...` | RefSeq assembly accession | NCBI RefSeq version of an assembly; do not additionally download if the same biological assembly has already been selected as a GCA |
| `SAMN_...` | NCBI BioSample | Links strain/sample metadata with NCBI records; not directly downloadable as an assembly FASTA |
| `SAMEA_...` | ENA BioSample | Corresponding BioSample type in the European archive |
| `PRJNA_...` | NCBI BioProject | Project container with many samples, reads, and/or assemblies |
| `PRJEB_...` | ENA BioProject | European project container |
| `SRR_...`, `ERR_...`, `DRR_...` | SRA/ENA/DDBJ run accessions | Raw reads; not primarily downloaded for the POC |
| `.fna` | FASTA genome file | Central input for Starfish after selecting suitable assemblies |
| `.gff`/`.gff3` | Genome annotation | Optional for gene/functional context and annotations |
| `.gbff` | GenBank flat file | Optional, comprehensive annotation container |

---

## 5. Technical Implementation

### 5.1 Computing Environment

The work is carried out under Windows with a Linux environment via WSL (Ubuntu). The bioinformatics environment resides in the Linux filesystem and not in the Windows-mounted area under `/mnt/c/...`.

Working directory:

```bash
~/promotion/barragan/data
```

This location is advantageous for bioinformatics analyses, since many small file operations and later mapping/assembly steps run considerably faster in the WSL filesystem than in the mounted Windows filesystem.

### 5.2 Software Environment

Installed and functional:

- Miniforge as a lightweight conda distribution
- Conda environment: `starfish_env`
- Python 3.8 within this environment
- Starfish as a tool for annotating large mobile elements and starship-like regions
- NCBI Datasets CLI: `datasets`
- NCBI metadata formatter: `dataformat`

Activating the working environment:

```bash
conda activate starfish_env
cd ~/promotion/barragan/data
```

The successful technical test was carried out via:

```bash
starfish --help
```

### 5.3 Exporting Assembly Metadata as TSV

The metadata fields used comprise:

| Field | Meaning in the POC |
|---|---|
| `Assembly Accession` | Download ID: GCA/GCF |
| `Organism Name` | Taxonomic designation; may read *Pyricularia oryzae* instead of *Magnaporthe oryzae* |
| `Assembly Name` | Submitter/assembly name; may contain a strain name |
| `Assembly Level` | Coarse quality tier: Contig, Scaffold, Chromosome, or Complete Genome |
| `Assembly Method` | Assembler used, e.g. Canu, HGAP, NextDenovo, Celera, or SPAdes |
| `Assembly Sequencing Tech` | Sequencing technology, e.g. Illumina, PacBio, or Oxford Nanopore |
| `Assembly BioProject Accession` | Project linkage |
| `Assembly BioSample Accession` | Central key for unambiguous merging with `TableS5` |
| `Assembly Stats Total Sequence Length` | Plausibility check for genome size |
| `Assembly Stats Number of Contigs` | Fragmentation measure; lower is better |
| `Assembly Stats Contig N50` | Contiguity measure; higher is better |
| `Assembly Stats Scaffold N50` | Additional contiguity measure |
| `Assembly Stats Genome Coverage` | Sequencing depth, where available |

---

## 6. Sequencing Technologies in the NCBI Dataset

A frequency analysis of the technology column showed that the overall dataset consists predominantly of Illumina assemblies, but includes a substantial proportion of long-read/hybrid assemblies.

Examples of the most frequently reported technologies:

| Technology | Number of records |
|---|---:|
| Illumina HiSeq | 134 |
| Illumina HiSeq 2000 | 116 |
| Illumina NovaSeq | 63 |
| Illumina HiSeq 2500 | 60 |
| Illumina MiSeq | 52 |
| Illumina GAIIx | 36 |
| Oxford Nanopore MinION | 27 |
| PacBio Sequel | 22 |
| Oxford Nanopore + Illumina | 12 |
| Oxford Nanopore (other variants) | at least 6 |
| PacBio (other variants) | at least 5 |
| PacBio/Illumina hybrid variants | several further records |

The technical labels are not fully standardized. For example, `minION`, `MiniION`, `Oxford Nanopore MinION`, `PromethION`, `PacBio RSII`, `PacBio Sequel`, `PacBio; Illumina MiSeq`, and further variants occur. For this reason, filtering was not performed on an exact term but with a broad, case-insensitive pattern.

Filter used:

```bash
grep -Ei 'pacbio|nanopore|minion|promethion|hybrid' \
  ncbi_m_oryzae_assemblies.tsv \
  > ncbi_m_oryzae_longread_candidates.tsv
```

The header line was then added:

```bash
{
  head -n 1 ncbi_m_oryzae_assemblies.tsv
  cat ncbi_m_oryzae_longread_candidates.tsv
} > ncbi_m_oryzae_longread_candidates_with_header.tsv
```

Result:

- **95 rows** in `ncbi_m_oryzae_longread_candidates_with_header.tsv`
- of which **1 header + 94 long-read/hybrid candidates**

---

## 7. Quality Criteria for Starfish Candidates

### 7.1 Rationale

Starships are large mobile elements, often in the range of several tens to several hundred kilobases. Accessory mini-chromosomes can likewise span large, repeat-rich, and structurally dynamic sequence regions. Heavily fragmented Illumina assemblies with thousands of contigs and N50 values in the kilobase range are therefore unsuitable for reliably determining element boundaries, structural variation, and integration context.

An observed example from the metadata set are Illumina GAIIx/Velvet assemblies with roughly 5,600–9,450 contigs and contig N50 values of about 7.6–16 kb. These can still be useful for global mapping or PAV questions, but do not form a sensible primary basis for structural Starship discovery.

### 7.2 Initial, Deliberately Inclusive Quality Filter

A long-read/hybrid assembly was included in the first high-quality candidate pool if it met at least one of the following criteria:

- `Assembly Level` is `Chromosome` or `Complete Genome`.
- Number of contigs is less than or equal to 200.
- Contig N50 is at least 1,000,000 bp.

The filter reads:

```bash
awk -F '\t' '
BEGIN { OFS="\t" }
NR == 1 { print; next }
$4 ~ /Chromosome|Complete Genome/ || ($12 != "" && $12 <= 200) || ($13 != "" && $13 >= 1000000) {
    print
}
' ncbi_m_oryzae_longread_candidates_with_header.tsv \
> ncbi_m_oryzae_longread_highquality.tsv
```

The criteria are deliberately inclusive so that suitable long-read assemblies are not excluded prematurely. The final selection is then made through prioritization and checks for biological relevance and possible duplicates.

### 7.3 Observed High-Quality Examples

The initial visual inspection yielded numerous very high-quality candidates, for example:

| Assembly | Technology | Level | Contigs | Contig N50 | Initial priority |
|---|---|---|---:|---:|---|
| `GCA_003015475.2` | PacBio Sequel | Contig | 13 | 5.50 Mb | High |
| `GCA_003015815.2` | PacBio Sequel | Contig | 13 | 5.45 Mb | High |
| `GCA_003015975.2` | PacBio Sequel | Contig | 10 | 5.64 Mb | High |
| `GCA_003016745.2` | PacBio Sequel | Contig | 16 | 6.16 Mb | High |
| `GCA_004346965.1` | PacBio RSII + Illumina | Complete Genome | 7 | 6.13 Mb | High |
| `GCA_004785725.2` | Oxford Nanopore | Chromosome | 10 | 6.47 Mb | High |
| `GCA_012272995.1` | Oxford Nanopore | Complete Genome | 9 | 6.16 Mb | High |
| `GCA_021442365.1` | PacBio + Illumina | Contig | 21 | 5.38 Mb | High |
| `GCA_021764705.1` | Oxford Nanopore + Illumina | Chromosome | 14 | 5.53 Mb | High |

These examples show that a structural reference panel for the Starfish POC is technically feasible.

### 7.4 Quality Classes for Later Selection

| Category | Criteria | Use |
|---|---|---|
| `high` | Long read/hybrid and Chromosome/Complete Genome, or at most 30 contigs, or contig N50 at least 5 Mb | Primary Starfish analysis; preferred references |
| `medium` | Long read/hybrid and at most 200 contigs or contig N50 at least 1 Mb | Panel supplementation; used after visual inspection |
| `exclude` | No clear long-read/hybrid basis or insufficient contiguity | Not used for primary structural detection; possibly usable later for mapping/PAV |

---

## 8. Special Quality Controls

### 8.1 Avoiding GCA/GCF Duplicates

A GCA and a GCF accession can represent the same underlying biological assembly. They must not be counted as independent genomes or downloaded twice for Starfish.

### 8.2 Identifying Repeated Submissions

Several assemblies with identical characteristics are visible in the initial candidate pool. Examples:

| Preferred assembly | Potentially redundant assembly | Evidence |
|---|---|---|
| `GCA_003015475.2` | `GCA_011799965.1` | same genome length, 13 contigs, contig N50 5.50 Mb |
| `GCA_003015815.2` | `GCA_011799905.1` | same genome length, 13 contigs, contig N50 5.45 Mb |
| `GCA_003015975.2` | `GCA_011799925.1` | same genome length, 10 contigs, contig N50 5.64 Mb |
| `GCA_003016745.2` | `GCA_011799915.1` | same genome length, 16 contigs, contig N50 6.16 Mb |

These candidates must be checked before the final selection via BioSample, assembly name, submitter information, and, if necessary, NCBI assembly details. Such redundancies must not be allowed to inflate the number of biologically independent isolates.

### 8.3 Accounting for Taxonomic Synonyms

The NCBI metadata may show `Pyricularia oryzae`, even though the PhD project uses the name *Magnaporthe oryzae*. This is a taxonomic synonym or a different nomenclatural convention and must not automatically be treated as grounds for exclusion during the data merge.

### 8.4 Not Blindly Trusting the Long-Read Label

A long-read technology entry is a strong selection criterion, but not complete proof of quality. Conversely, a high-quality assembly may exist despite an empty technology field. Prioritization is therefore always carried out in combination via:

- sequencing technology
- assembly level
- contig count
- contig N50
- scaffold N50
- genome length and coverage as plausibility values
- biological metadata and linkage with the Barragán dataset

---

## 9. Current Data Holdings and Folder Structure

Recommended/used project structure:

```text
~/promotion/barragan/
├── data/
│   ├── Copy of Barragan2024_SupplementalTables.xlsx
│   ├── ncbi_m_oryzae_assemblies.jsonl
│   ├── ncbi_m_oryzae_assemblies.tsv
│   ├── ncbi_m_oryzae_longread_candidates.tsv
│   ├── ncbi_m_oryzae_longread_candidates_with_header.tsv
│   └── ncbi_m_oryzae_longread_highquality.tsv
├── metadata/
│   ├── barragan_longread_master.xlsx
│   ├── assembly_accessions_selected.txt
│   └── download_manifest.tsv
├── assemblies/
│   ├── ncbi_dataset/
│   ├── genome_fasta/
│   ├── annotations/
│   └── sequence_reports/
├── starfish/
│   ├── input/
│   ├── output/
│   └── logs/
└── scripts/
    ├── build_master_metadata.py
    ├── select_assemblies.py
    └── run_starfish.sh
```

To date, only metadata files have been generated. No FASTA, GFF3, protein, or FASTQ files have been collected.

---

## 10. Next Steps

### Step 1: Generate the Master Spreadsheet

The next step planned immediately is to link:

- `TableS5` from `Copy of Barragan2024_SupplementalTables.xlsx`, and
- `ncbi_m_oryzae_longread_highquality.tsv`.

The primary join key is:

```text
TableS5: Sample  <->  NCBI: Assembly BioSample Accession
```

The master file should retain all Barragán isolates and document the match status for each isolate:

| Status | Meaning |
|---|---|
| `matched` | Matching high-quality long-read/hybrid assembly found via BioSample |
| `no_highquality_longread_match` | No match in the high-quality long-read pool; does not necessarily mean that no raw data or assembly exists |
| `manual_check` | Uncertain or ambiguous match; manual review required |

Target filename:

```text
barragan_longread_master.xlsx
```

### Step 2: Select the Final Starfish Panel

An initial structural panel of roughly 20–40 biologically independent assemblies will be compiled from the master file. Selection criteria:

- Representation of different lineages: Oryza, Triticum, *Eleusine*, Lolium, and other wild-grass lineages.
- Highest quality takes priority over largest number.
- No GCA/GCF duplicates and no multiply submitted identical assemblies.
- As many different geographic and host-related contexts as possible.
- Inclusion of isolates related to Barragán studies, in particular known mChrA/mChr carriers, where the data are publicly available.
- Additionally, reference isolates such as 70-15, provided a high-quality matching assembly is available.

### Step 3: Download the Finally Selected Assemblies

Only after the selection is a file `assembly_accessions_selected.txt` generated, with exactly one `GCA_...` or `GCF_...` accession per line.

Example:

```text
GCA_003015475.2
GCA_004346965.1
GCA_004785725.2
```

FASTA, GFF3, proteins, GenBank flat files, and sequence reports are then downloaded:

```bash
./datasets download genome accession \
  --inputfile assembly_accessions_selected.txt \
  --include genome,gff3,protein,gbff,seq-report \
  --filename m_oryzae_starfish_panel.zip

mkdir -p ../assemblies/ncbi_dataset
unzip m_oryzae_starfish_panel.zip -d ../assemblies/ncbi_dataset
```

The genomic FASTA files can then be checked with:

```bash
find ../assemblies/ncbi_dataset -type f -name "*.fna" | sort
```

### Step 4: Standardize FASTA Files and Quality Control

Before Starfish:

- Unambiguous assignment of each FASTA to GCA/GCF, strain ID, BioSample, and lineage.
- Checking FASTA file size, contig count, and headers.
- Documentation of exclusions.
- Separate storage of genome FASTA, annotation, and metadata.

### Step 5: Starfish Test Run on 1–3 Reference Assemblies

No full analysis across the entire panel initially. A controlled test run is recommended on:

- one PacBio assembly with very high contiguity,
- one ONT/hybrid assembly,
- one biologically relevant reference or a known mChr-associated isolate.

This will check:

- whether Starfish runs technically reproducibly;
- whether Captain candidates and starship-like regions are output;
- whether the results are plausible in terms of length, repeat context, and contig position;
- which output formats are produced for the subsequent comparison and network phase.

### Step 6: Starship Catalog and Sequence Comparison

After a successful test, a candidate catalog will be generated, at minimum with:

| Field | Content |
|---|---|
| `element_id` | Unique candidate identifier |
| `assembly_accession` | Source assembly |
| `strain_id` | Strain, where resolvable |
| `lineage` | Host/lineage assignment |
| `element_type` | Starship, starship-like, mChr-associated, unclear |
| `contig` | Carrier contig |
| `start`, `end` | Coordinates |
| `length_bp` | Length |
| `captain_gene` | Captain/tyrosine recombinase evidence |
| `boundary_evidence` | e.g. TSD, terminal sequences, flanking structure |
| `quality_class` | High/Medium/Low |

In the next step, the elements will be compared pairwise across lineages. For the POC, a pragmatic criterion is suitable:

- at least 95% nucleotide identity
- over at least 80% of the element length
- occurrence in at least two different host-associated lineages

These cases are not yet automatically demonstrated as HGT. They are **prioritized HGT candidates**, which must subsequently be checked against the core-genome phylogeny, regional similarity, PAV, and, where applicable, D-statistics.

### Step 7: Transition to the Global PAV/Mapping Level

Only the structurally well-defined candidates from the high-quality assemblies subsequently serve as references for the large global short-read dataset. This makes the following question answerable:

> In which lineages, hosts, countries, and isolates are the same or highly similar mobile elements present?

This is the point at which the coverage/mChr information already present in the Barragán dataset (TableS6) and the global BioSample/BioProject linkages become particularly relevant.

---

## 11. Decision Criteria for the POC

| POC question | Positive outcome | Consequence |
|---|---|---|
| Are there enough high-quality assemblies? | At least approximately 20 biologically independent long-read/hybrid assemblies across several lineages | Build the Starfish panel and start the structural analysis |
| Are there starship-like/mChr-associated candidates? | Reproducible Starfish candidates with Captain/structural features | Build the candidate catalog and sequence comparison |
| Are there highly similar elements across lineages? | Elements with high identity and high coverage in several lineages | Move candidates into the HGT evidence pipeline |
| Are the events sufficiently frequent? | Several independent donor–recipient/element relationships | Prepare local recovery approximation and GLMM power analysis |
| Are the events too rare? | Only isolated cases or no recurring relationships | Simplify the model; pool mChr/Starships; prioritize exploratory goals; consider targeted new sampling |

---

## 12. Conclusion on the Current Status

The workflow to date answers the first technical feasibility question positively:

- The NCBI query yields 605 public assembly records for *Magnaporthe/Pyricularia oryzae*.
- Among these, 94 long-read or hybrid candidates were identified based on sequencing technology.
- The quality filter reveals numerous structurally very suitable assemblies with 7–56 contigs and contig N50 values of approximately 3–6.5 Mb.
- A qualitatively robust, assembly-based Starfish reference panel is therefore realistic for a POC.
- The next critical work step is not the full download, but the controlled linkage of these assemblies with the biological metadata from Barragán TableS5 and the removal of biological or technical duplicates.

The work carried out so far should thus be classified as a reproducible data-foundation and quality phase. HGT has not yet been demonstrated, and the observed assembly availability must not yet be interpreted as a frequency of horizontal transfers.
