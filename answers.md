---
title: "Reproducible Sequence Analysis Pipeline"
subtitle: "HPDC: Assignment 3"
author: "David Hidalgo Fàbregas"
date: "15 june 2026"
geometry: margin=2.5cm
lang: en-US
---


# 1. Dynamic Programming Edit Distance

### Complete Dynamic Programming matrix for the first sequence in the FASTA file

### What each cell $M(i,j)$ represents?
Each cell $M(i,j)$ represents the minimum edit distance required to transform the first $i$ characters of the reference sequence (pattern[0:i]) into the first $j$ characters of the target sequence (text[0:j]). Finally, it stores the optimal solution to each smaller subproblem.

### Why Dynamic Programming avoids the repeated computations that appear in a brute-force recursive solution?
In a recursive brute-force solution (top-down without memory), the algorithm recalculates the distance between the same sequence fragments over and over again, resulting in exponential runtime. Dynamic Programming, using tabulation (bottom-up), solves each subproblem only once and stores the result in the matrix. When it needs that value to calculate larger alignments, it simply retrieves it from memory, avoiding redundant work.

### Why the time complexity of this tabulation algorithm is $O(n \cdot m)$?
The algorithm consists of filling in a matrix of dimensions $(n+1) \times (m+1)$. To calculate the value of each individual cell, the algorithm performs a constant number of operations, $O(1)$: checking for character matches and finding the minimum among three adjacent values (diagonally, above, and to the left). Since $O(1)$ operations are performed for each of the $n \cdot m$ cells, the total time complexity is strictly proportional to the area of the matrix, that is, $O(n\cdot m)$.

# 2. Regular Expressions Analysis of CIGAR Strings

### Why Regular Expressions are useful for extracting information from CIGAR strings and other bioinformatis file formats?
Bioinformatics formats such as CIGAR strings condense complex information about sequence alignments into text strings by combining letters that represent operations. Regular expressions provide a declarative, fast, and very powerful language for searching for, isolating, and counting specific patterns within these strings without the need to program complex manual iterative loops.

# 3. Conda Environment

### Which command creates the environment?

```bash
conda create -n hpdc_env python=3 snakemake -c bioconda -c conda-forge -y
```
The `-c` flag defines the channels where the packages can be found.

The `-y` flag automatically answers “yes” to any confirmation prompts that appear during the installation process.

### Which command activates the environment?
```bash
conda activate hpdc_env
```

### Which command exports the environment?
```bash
conda env export > environment.yml
```
With only one `>` instead of two `>>`, the command overwrites an existing file, if it exists (if not, it creates a new file).

This file will only be compatible with machines running the same operating system as the original machine (Linux) because it also exports the binary compilation code, which depends on the processor architecture and operating system. Appart from that, Snakemake is only available for Linux/macOS.

### Why is an environment file useful for reproducibility?
An environment file is essential for reproducibility because it guarantees that an analysis can be executed on different computers and at different times with the exact same software configuration. Without it, using different versions of Python or software packages could lead to inconsistent results or executions failures. By tracking the packages versions used and sharing this file, any user can recreate the isolated software environment perfectly, ensuring the workflow runs exactly as intended.

### Which packages were included in your environment?
Only the python=3 and snakemake packages were included.

### Why these packages are required for the workflow?

- **Python 3** is the base interpreter required to run all the analysis scripts (parse_fasta.py, alignment.py, cigar_stats.py), which includes `sys` and `re` (regex) libraries.

- **Snakemake** is the workflow manager required to link scripts, manage their dependencies, and automate the execution of the pipeline in a reproducible manner.

### How another user could recreate your environment usin the exported `environment.yml` file?
Any other user who receives the environment.yml file simply needs to have Conda installed and run the following command in their terminal:

```bash
conda env create -f environment.yml
```

This will read the exact dependencies and versions and automatically build an identical environment.


# 4. Workflow Automation with Snakemake

### Output of `snakemake -n`

### Generated workflow DAG.

### What each rule does?

### Why Snakemake improves reproducibility compared with manually running each script?