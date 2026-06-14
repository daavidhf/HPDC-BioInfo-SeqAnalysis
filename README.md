# Reproducible Sequence Analysis Pipeline (Assignment 3)

This repository contains the complete workflow for **Assignment 3**, which focuses on building a reproducible bioinformatics pipeline from scratch. The project demonstrates core skills in algorithm implementation, regular expressions, environment management, and workflow automation.

## 📝 Assignment Overview

The main objective of this assignment was to process a set of biological sequences provided in a FASTA file, perform local alignments using Dynamic Programming, analyze the alignment strings (CIGAR), and ensure the entire process is 100% reproducible using **Conda** and **Snakemake**.

## 🚀 Steps Followed (Tasks)

The project was developed following a strict Git branching strategy, completing the following tasks sequentially:

* **Task 1: FASTA Parsing (`parse_fasta.py`)**
    * Developed a script to read `PM_50.fasta`, extract sequence accessions, and calculate the length of each sequence.
    * Output: `sequences.tsv`

* **Task 2: Dynamic Programming Alignment (`alignment.py`)**
    * Implemented a bottom-up Dynamic Programming algorithm ($O(n \cdot m)$ complexity) to compute the minimum edit distance.
    * Aligned the first 18 nucleotides of each sequence against the first 18 nucleotides of the first sequence (used as the reference pattern).
    * Output: `distances.tsv`

* **Task 3: CIGAR String Generation (Traceback)**
    * Extended `alignment.py` to include a traceback algorithm.
    * Walked backward through the DP matrix to generate expanded CIGAR strings (Matches `M`, Substitutions `X`, Insertions `I`, Deletions `D`).
    * Output: `alignments.tsv`

* **Task 4: Regex Analysis and Reporting (`cigar_stats.py`)**
    * Used Python's `re` module (Regular Expressions) to parse the CIGAR strings and count the occurrences of each operation.
    * Extracted the sequences with the minimum/maximum edit distance and the highest number of matches.
    * Outputs: `cigar_stats.tsv` and `report.txt`

* **Task 5: Reproducible Environment (`environment.yml`)**
    * Created an isolated Conda environment containing only the strictly necessary dependencies (`python=3` and `snakemake`).
    * Exported the configuration using `--from-history` for cross-platform compatibility.

* **Task 6: Workflow Automation (`Snakefile`)**
    * Wrote a Snakemake workflow to connect all the independent Python scripts into a Directed Acyclic Graph (DAG).
    * Automated the execution order based on input/output dependencies, eliminating the need to run scripts manually.

---

## 📁 Repository Structure

```text
assignment3/
├── .gitignore             # Excludes Python caches and Snakemake metadata
├── README.md              # Project documentation
├── Snakefile              # Snakemake workflow configuration
├── environment.yml        # Conda environment definition
├── PM_50.fasta            # Input raw biological sequences
├── parse_fasta.py         # Script for Task 1
├── alignment.py           # Script for Tasks 2 & 3
├── cigar_stats.py         # Script for Task 4
├── dag.png                # Visual representation of the workflow DAG
├── answers.pdf            # Theoretical questions and answers for the assignment
└── [Outputs]              # Generated files (sequences.tsv, distances.tsv, etc.)