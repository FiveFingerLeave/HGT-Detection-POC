# Starfish Setup and First Successful YR Test

**Project:** PhD thesis – Barragan / *Magnaporthe oryzae*  
**Date:** August 31, 2026  
**Purpose of this document:** Reproducible documentation of the steps carried out so far to prepare and first test Starfish.

> This document deliberately distinguishes between successfully executed steps, warnings, and work that is still open. It documents the current status; it does not yet claim that Starship elements or mini-chromosomes have been conclusively identified.

---

## 1. Scientific Objective

This PhD project investigates genomic variation, horizontal gene transfer (HGT), starship-like mobile elements, and accessory mini-chromosomes in isolates of *Magnaporthe oryzae*.

Starfish is a tool for detecting and annotating large mobile elements, in particular so-called **Starships**. A central feature of many Starships is a gene encoding a tyrosine recombinase (YR), often referred to as the **Captain gene**. A single YR hit, however, is initially only a candidate and not complete proof of a Starship element.

The analysis carried out so far therefore had the limited goal of:

1. Testing Starfish in the existing conda environment.
2. Preparing an NCBI assembly as input.
3. Searching for YR candidates in this assembly via de novo gene prediction and HMM validation.

---

## 2. Working Environment and Directories

Conda environment used:

```bash
starfish_env
```

Key working directories:

```text
~/promotion/barragan/starfish
~/promotion/barragan/data
~/promotion/barragan/assemblies/ncbi_dataset/ncbi_dataset/data
```

Test assembly used:

```text
GCA_004346965.1
Assembly name: ASM434696v1
Genome FASTA:
~/promotion/barragan/assemblies/ncbi_dataset/ncbi_dataset/data/GCA_004346965.1/GCA_004346965.1_ASM434696v1_genomic.fna
```

---

## 3. Generating the NCBI Assembly Report

First, a tab-separated table was created from the NCBI JSONL assembly report:

```bash
cd ~/promotion/barragan/data

./dataformat tsv genome \
  --inputfile "$HOME/promotion/barragan/assemblies/ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl" \
  --fields accession,assminfo-biosample-accession,assminfo-name \
  > "$HOME/promotion/barragan/starfish/input/assembly_data_report.tsv"
```

The first rows looked like:

```text
Assembly Accession	Assembly BioSample Accession	Assembly Name
GCA_000292605.2	SAMN02981399	PoP131
GCA_003016745.2	SAMN06050113	ASM301674v2
GCA_004346965.1	SAMN10491321	ASM434696v1
```

### Result and Problem

The attempt to use

```bash
starfish format-ncbi \
  --report input/assembly_data_report.tsv \
  --assemblies input/test_assemblies.txt
```

failed. Starfish could not even parse the header row as a data row.

A headerless file with two columns was then created:

```text
GCA_004346965.1	ASM434696v1
```

The tabs were confirmed with the following test:

```bash
grep '^GCA_004346965.1' input/assembly_data_report_starfish.tsv | cat -A
```

Output:

```text
GCA_004346965.1^IASM434696v1$
```

`^I` means: an actual tab character.

Despite correct tabs, `starfish format-ncbi` rejected the file. From this it follows that:

- The problem was not tab-separation.
- `format-ncbi` presumably expects a more specific report structure than the two-column table generated.
- This subcommand was not pursued further for the first functional Starfish test.

The file can be kept, but is not required for the direct `annotate` approach documented below.

---

## 4. Preparing the Assembly File for Starfish

The help output of `starfish annotate` shows that `--assembly` expects a two-column TSV file:

```text
genomeID<TAB>path-to-assembly-FASTA
```

The following file was therefore created:

```bash
cd ~/promotion/barragan/starfish

printf '%s\t%s\n' \
"GCA_004346965.1" \
"$HOME/promotion/barragan/assemblies/ncbi_dataset/ncbi_dataset/data/GCA_004346965.1/GCA_004346965.1_ASM434696v1_genomic.fna" \
> input/test_assemblies_2col.tsv
```

Contents, checked with visible control characters:

```bash
cat input/test_assemblies_2col.tsv | cat -A
```

Output:

```text
GCA_004346965.1^I/home/flori/promotion/barragan/assemblies/ncbi_dataset/ncbi_dataset/data/GCA_004346965.1/GCA_004346965.1_ASM434696v1_genomic.fna$
```

This file is correctly formatted: the first column contains the stable genome ID, the second the full FASTA path.

---

## 5. Existing Starfish Databases

The following database files were found in the conda environment:

```bash
find "$CONDA_PREFIX/db" -type f | sort
```

```text
YRsuperfamRefs.faa
YRsuperfams.p1-512.hmm
duf3723.hmm
duf3723.mycoDB.faa
fre.hmm
fre.mycoDB.faa
myb.SRG.fa
myb.hmm
nlr.hmm
nlr.mycoDB.faa
plp.hmm
plp.mycoDB.faa
```

For the first Starship-oriented search run, the following were used:

```text
YRsuperfamRefs.faa          Reference amino acid sequences for YR genes
YRsuperfams.p1-512.hmm      HMM profiles for validating YR superfamilies
```

An **HMM** (Hidden Markov Model) is a statistical profile of a protein family. It recognizes conserved sequence patterns even when the proteins are not exactly identical. Combining reference protein search with HMM validation reduces nonspecific hits.

---

## 6. Local Help for `starfish annotate`

The local help confirmed the key required arguments:

```text
-a, --assembly    2 column tsv: genomeID, path to assembly FASTA
-p, --profile     profile HMM file
-P, --proteins    FASTA file of query amino acid sequences
-x, --prefix      prefix for naming all output files
-i, --idtag       prefix for predicted gene featureIDs
-o, --outdir      output directory
```

Optionally, `-g` can be used to pass a two-column table of genome ID and GFF file path. This will later be important for merging existing NCBI gene models with newly predicted YR genes.

---

## 7. Successful First YR Test Run

Before the run, output, temporary, and log directories were prepared:

```bash
mkdir -p output/test_YR temp/test_YR logs
```

Command executed:

```bash
starfish annotate \
  --assembly input/test_assemblies_2col.tsv \
  --profile "$CONDA_PREFIX/db/YRsuperfams.p1-512.hmm" \
  --proteins "$CONDA_PREFIX/db/YRsuperfamRefs.faa" \
  --prefix GCA_004346965.1_YR \
  --idtag YR \
  --outdir output/test_YR \
  --tempdir temp/test_YR \
  --threads 6 \
  > logs/annotate_GCA_004346965.1_YR.log 2>&1
```

### Result

The run completed successfully. Relevant log lines:

```text
running metaeuk easy-predict for 1 assemblies..
running hmmsearch on metaeuk annotations..
filtering metaeuk annotations based on hmmsearch results..
found 13 new YR genes and 0 YR genes that overlap with 0 existing genes
done
```

**Interpretation:** In the assembly `GCA_004346965.1`, 13 newly predicted genes were found whose proteins satisfied the YR HMM criterion used. These genes are **YR candidates**. They are possible Captain genes, but not yet confirmed Starship elements.

### Generated Result Files

```text
output/test_YR/GCA_004346965.1_YR.fas
output/test_YR/GCA_004346965.1_YR.filt.fas
output/test_YR/GCA_004346965.1_YR.filt.gff
output/test_YR/GCA_004346965.1_YR.filt.ids
output/test_YR/GCA_004346965.1_YR.filt.old2new.ids
output/test_YR/GCA_004346965.1_YR.gff
output/test_YR/GCA_004346965.1_YR.hmmout
output/test_YR/GCA_004346965.1_YR.hmmout.ids
output/test_YR/GCA_004346965.1_YR.metaeuk.log
```

### Important Files

| File | Meaning |
|---|---|
| `*.gff` | Genomic coordinates of the MetaEuk-predicted genes |
| `*.filt.gff` | YR gene models filtered after HMM validation; more important than the unfiltered GFF for the next steps |
| `*.fas` | Amino acid or sequence output of the predicted genes |
| `*.filt.fas` | Sequences of the filtered, HMM-validated YR hits |
| `*.hmmout` | Detailed output of `hmmsearch` |
| `*.hmmout.ids` | IDs of the HMM hits |
| `*.filt.ids` | IDs of the final accepted YR candidates |
| `*.metaeuk.log` | Detailed log of the de novo gene prediction |

---

## 8. Header Warnings

During the successful run, warnings appeared such as:

```text
warning: CP034210.1 ... is being parsed into <2 components using separator '_'.
Make sure ALL sequence headers are formatted like <genomeID><separator><featureID>
```

The original FASTA headers apparently contain only the contig or chromosome ID:

```text
>CP034210.1
```

For multi-genome workflows, however, Starfish expects a unique composite identifier, for example:

```text
>GCA_004346965.1_CP034210.1
```

### Why This Matters

The same contig names can occur across several isolates. If, for example, a contig `contig_1` or `CP034210.1` exists in each assembly, the origin would no longer be unambiguous later on. Prefixing with the genome ID makes each sequence globally unique.

### Consequence

The test run is valid as proof that Starfish works. For the systematic multi-isolate pipeline, however:

1. FASTA headers should be standardized.
2. Associated GFF files should be adjusted to use exactly the same sequence IDs.
3. The YR annotation should then be run again.

---

## 9. What Has Been Established So Far

- Starfish can be run in the conda environment `starfish_env`.
- Direct submission of a two-column assembly TSV to `starfish annotate` works.
- MetaEuk and HMMER ran successfully on the test assembly.
- 13 HMM-validated YR candidates were found in `GCA_004346965.1`.
- The original NCBI FASTA headers should be normalized for the multi-isolate analysis.
- `starfish format-ncbi` has not yet been successfully set up and is not required for the current direct analysis path.

---

## 10. What Has Not Yet Been Established

The following statements **must not** be inferred from the test conducted so far:

- 13 complete Starship elements have been discovered.
- The 13 YR genes actually lie within mobile elements.
- The candidates possess cargo genes or well-defined element boundaries.
- A candidate is located on a mini-chromosome.
- Horizontal transfers have been demonstrated between isolates.

These statements require the following workflow steps: standardization of all inputs, consolidation with complete gene annotations, Starship context analysis, mini-chromosome classification, and cross-isolate comparisons.

---

## 11. Reproducibility

For every future run, at least the following information should be stored:

- Exact command executed
- Tool and conda environment versions
- Input file and checksum
- Parameters, in particular HMM e-value, threads, and MetaEuk options
- Complete logs
- Date, analysis purpose, and git commit of the workflow

The log file of the successful first test is located at:

```text
logs/annotate_GCA_004346965.1_YR.log
```

The current results are located at:

```text
output/test_YR/
```

---

## 12. Next Documented Work Step

Next, a structured, automated pipeline will be created. It should:

1. harmonize all FASTA and GFF IDs,
2. run Starfish reproducibly per isolate,
3. transparently classify Starship candidates,
4. independently identify and compare mini-chromosomes,
5. evaluate presence/absence, sequence similarity, and phylogenetic discordance between isolates.

A detailed guide for this is provided in the separate document **"Starfish_and_MiniChromosome_Analysis_Plan.md"**.
