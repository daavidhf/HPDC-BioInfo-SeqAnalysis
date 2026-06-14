# =====================================================
# Assignment 3: Reproducible Sequence Analysis Pipeline
# =====================================================


# Target rule: Defines the final file we want to generate.
rule all:
    input:
        "report.txt",
        "sequences.tsv"

# Part 1: Extract sequences and compute lengths from the FASTA file
rule parse:
    input:
        "PM_50.fasta"
    output:
        "sequences.tsv"
    shell:
        "python parse_fasta.py {input} {output}"

# Parts 2 & 3: Compute Edit Distance and generate CIGAR strings using DP
rule align:
    input:
        "PM_50.fasta"
    output:
        dist="distances.tsv",
        aln="alignments.tsv"
    shell:
        "python alignment.py {input} {output.dist} {output.aln}"

# Part 4: Analyze CIGAR strings using Regular Expressions and generate summary
rule stats:
    input:
        "alignments.tsv"
    output:
        stats="cigar_stats.tsv",
        report="report.txt"
    shell:
        "python cigar_stats.py {input} {output.stats} {output.report}"
