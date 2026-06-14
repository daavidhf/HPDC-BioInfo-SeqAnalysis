import sys
import re

def edit_distance_dp(pattern, text):
    """
    Computes the edit distance between two sequences using Dynamic Programming.
    pattern: The reference sequence fragment (length n).
    text: The target sequence fragment being compared (length m).
    """
    n = len(pattern)
    m = len(text)

    # Initialize the DP matrix with zeros. Size: (n+1) x (m+1)
    # M(i, j) represents the edit distance between pattern[0:i] and text[0:j]
    dp_matrix = [[0] * (m + 1) for _ in range(n + 1)]

    # Initialize the first column (transforming pattern of length i to an empty string requires i deletions)
    for i in range(n + 1):
        dp_matrix[i][0] = i  # Deletion cost
    
    # Initialize the first row (transforming an empty string to text of length j requires j insertions)
    for j in range(m + 1):
        dp_matrix[0][j] = j  # Insertion cost

    # Fill the DP matrix using bottom-up tabulation
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Check if nucleotides match (cost 0) or substitute (cost 1)
            if pattern[i-1] == text[j-1]:
                cost = 0
            else:
                cost = 1
            
            # Calculate minimum cost among the three possible movements:
            diagonal = dp_matrix[i-1][j-1] + cost   # Match/Substitution
            up = dp_matrix[i-1][j] + 1              # Deletion
            left = dp_matrix[i][j-1] + 1            # Insertion
            dp_matrix[i][j] = min(diagonal, up, left)
    
    # Final edit distance is located at the bottom-right cell of the matrix
    final_distance = dp_matrix[n][m]
    return dp_matrix, final_distance

def get_cigar(pattern, text, dp_matrix):
    """
    Performs a traceback on the DP matrix to generate the expanded CIGAR string.
    """
    i = len(pattern)
    j = len(text)
    cigar_ops = []

    # Traceback from bottom-right to top-left
    while i > 0 or j > 0:
        # Check diagonal movement first (match of substitution)
        if i > 0 and j > 0:
            cost = 0 if pattern[i-1] == text[j-1] else 1
            if dp_matrix[i][j] == dp_matrix[i-1][j-1] + cost:
                if cost == 0:
                    cigar_ops.append('M')  # Match
                else:
                    cigar_ops.append('S')  # Substitution
                i -= 1
                j -= 1
                continue
        
        # Check upward movement (deletion)
        if i > 0 and dp_matrix[i][j] == dp_matrix[i-1][j] + 1:
            cigar_ops.append('D')  # Deletion
            i -= 1
            continue

        # Check left movement (insertion)
        if j > 0 and dp_matrix[i][j] == dp_matrix[i][j-1] + 1:
            cigar_ops.append('I')  # Insertion
            j -= 1
            continue

    return "".join(reversed(cigar_ops))

def print_dp_matrix(pattern, text, matrix):
    """
    Prints the complete Dynamic Programming matrix in a readable format to the console.
    This output can be copied directly into the answers.pdf report.
    """
    print("\n=============================================================")
    print("DYNAMIC PROGRAMMING MATRIX FOR THE FIRST SEQUENCE")
    print("=============================================================")

    # Print column headers (text characters)
    header = "      -  " + "  ".join(text)
    print(header)

    # Print each row of the matrix
    for i, row in enumerate(matrix):
        row_label = pattern[i-1] if i > 0 else "-"
        row_str = f" {row_label} | " + " ".join(f"{val:2d}" for val in row)
        print(row_str)
    print("=============================================================\n")

def parse_fasta_fragments(fasta_path):
    """
    Parses the FASTA file and reconstructs complete sequences.
    Returns a list of tuples containing (accession, full_sequence).
    """
    sequences = []
    current_accession = None
    current_seq_lines = []
    
    regex_accession = re.compile(r'^>\s*(\S+)')
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_accession:
                    sequences.append((current_accession, "".join(current_seq_lines)))
                
                match = regex_accession.search(line)
                current_accession = match.group(1) if match else "Unknown"
                current_seq_lines = []
            else:
                current_seq_lines.append(line)
                
        # Append the last sequence in the file
        if current_accession:
            sequences.append((current_accession, "".join(current_seq_lines)))
            
    return sequences

def main():
    # Setup command line arguments with default values required by the assignment
    input_fasta = sys.argv[1] if len(sys.argv) > 1 else "PM_50.fasta"
    out_distances = sys.argv[2] if len(sys.argv) > 2 else "distances.tsv"
    out_alignments = sys.argv[3] if len(sys.argv) > 3 else "alignments.tsv"
    
    try:
        # Step 1: Parse all sequences from the FASTA file
        sequences = parse_fasta_fragments(input_fasta)
        if not sequences:
            print("Error: No sequences found in the input file.")
            return
            
        # Step 2: Extract the reference pattern from the first sequence (first 18 bp)
        first_accession, first_full_seq = sequences[0]
        reference_pattern = first_full_seq[:18]
        
        # Step 3: Compute distances and write to the output TSV file
        with open(out_distances, 'w') as f_dist, open(out_alignments, 'w') as f_aln:
            f_dist.write("Accession\tDistance\n")
            f_aln.write("Accession\tDistance\tCIGAR\n")

            for idx, (accession, full_seq) in enumerate(sequences):
                # Extract the first 18 nucleotides of the current sequence
                target_text = full_seq[:18]
                
                # 1. Compute the DP matrix and the final minimum edit distance
                matrix, distance = edit_distance_dp(reference_pattern, target_text)
                # 2. Perform traceback to get CIGAR
                cigar_string = get_cigar(reference_pattern, target_text, matrix)
                
                # Write record to TSV
                f_dist.write(f"{accession}\t{distance}\n")
                f_aln.write(f"{accession}\t{distance}\t{cigar_string}\n")

                # Requirement: Print the complete DP matrix for the first sequence
                if idx == 0:
                    print_dp_matrix(reference_pattern, target_text, matrix)
                    print(f"First sequence CIGAR: {cigar_string}")
                    
        print(f"Success: Analysis completed.")
        print(f" - Distances saved to '{out_distances}'")
        print(f" - Alignments saved to '{out_alignments}'")
        
    except FileNotFoundError:
        print(f"Error: The file '{input_fasta}' was not found.")

if __name__ == "__main__":
    main()