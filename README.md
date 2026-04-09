# Mutation Counting Tool

A bioinformatics tool for analyzing mutational signatures from multiple sequence alignments (MSA). 
Given a FASTA alignment file, it computes a consensus sequence, counts trinucleotide (triplet) contexts, and categorizes point mutations using a strand-symmetric conversion matrix.

---

## Requirements

- Python 3.7+
- pandas
- tkinter (usually included with Python)
- `conversion_matrix.csv` (must be present in the working directory)
- `conversion_matrix.py` (must be in the same directory as `Mutation_counting.py`)

Install dependencies:

```bash
pip install pandas
```

---

## Files

| File | Description |
|---|---|
| `Mutation_counting.py` | Main script — runs the full analysis pipeline |
| `conversion_matrix.py` | Helper module for strand-symmetric mutation categorization |
| `conversion_matrix.csv` | Pre-computed 64×64 reverse-complement mapping matrix |

---

## Usage

Run the main script directly:

```bash
python Mutation_counting.py
```

A file dialog will open prompting you to select a `.fasta` file containing your multiple sequence alignment. All output files are saved to the same directory as the input file.

---

## Input Format

The input must be a **FASTA-formatted multiple sequence alignment (MSA)**:

- All sequences must be the **same length** (i.e., already aligned, with gap characters `-` where needed)
- Valid nucleotide characters: `A`, `T`, `C`, `G`, `U` (U is converted to T automatically)
- Any non-standard characters are replaced with `-`
- Both DNA and RNA sequences are accepted

Example:

```
>Virus_1
ATGCTAGCTAGCT---ATGC
>Virus_2
ATGCTAGCTAGCTATGATGC
```

---

## Outputs

All output files are written to the same directory as the input FASTA, using the input filename as a prefix.

| Output File | Description |
|---|---|
| `<name>_consensus_sequence.txt` | Consensus sequence computed from the MSA (threshold: 60% majority, minimum 5 sequences supporting a position) |
| `<name>_triplet_counts.csv` | Per-sequence counts of all trinucleotide contexts observed (rows = sequences, columns = triplets) |
| `<name>_mutation_counts.csv` | Per-sequence counts of mutations grouped into strand-symmetric categories (rows = sequences, columns = category IDs) |
| `<name>_mutation_category.csv` | Maps each category ID to the specific trinucleotide substitutions it represents |

---

## How It Works

### 1. Consensus Sequence
A position-wise consensus is computed from the MSA. A nucleotide is included in the consensus only if it appears in ≥60% of sequences at that position and is supported by at least 5 sequences. Positions that don't meet this threshold are marked as `-`.

### 2. Triplet Context Counting
For each sequence, triplet (trinucleotide) contexts are counted using the consensus sequence as a reference backbone. Gaps in a sequence are skipped; the consensus nucleotide fills those positions. This ensures triplet counts reflect the underlying genomic composition rather than alignment artifacts.

### 3. Mutation Categorization
Point mutations are identified by comparing each sequence to the consensus. For each mutation, the surrounding trinucleotide context (one base on each side) is extracted, and the mutation is assigned to a **strand-symmetric category** using the conversion matrix. This collapses reverse-complement equivalent mutations into a single category (e.g., C>T in context ACG and its reverse complement CGT>CAT are counted together), following standard mutational signature conventions.

---

## Notes

- The tool requires a proper MSA as input — raw unaligned sequences will trigger an error.
- The consensus is used as a proxy reference; sequences are not compared against an external reference genome.
- Positions at the very start and end of the alignment (index 0 and last) are excluded from mutation counting to avoid edge cases with triplet extraction.
- `conversion_matrix.csv` must be present in the **current working directory** when the script is run, not necessarily the same directory as the script itself.
