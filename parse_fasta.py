import sys
import re

def parse_fasta(input_file, output_file):
    # Regular expression to match FASTA headers
    # Search for lines starting with '>' (^>)
    # Ignore any optional whitespaces (\s*)
    # Capture first character block without spaces (\S+)
    regex_accession = re.compile(r'^>\s*(\S+)')
    current_accession = None
    current_length = 0

    try:
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            f_out.write("Accession\tLength\n")  # Write header to output file
            for line in f_in:
                line = line.strip() # Remove leading/trailing whitespace
                if line.startswith('>'):
                    # If we were reading a sequence, save the data
                    if current_accession is not None:
                        f_out.write(f"{current_accession}\t{current_length}\n")
                    
                    # Extract accession number using regex
                    match = regex_accession.search(line)
                    if match:
                        current_accession = match.group(1) # Extract the captured group (accession)
                    else:
                        current_accession = "Unknown"  # Default value if no accession is found

                    # Reset the length for the new sequence
                    current_length = 0
                else:
                    # If not a header line, it is DNA. Sum the length of the line.
                    current_length += len(line)

            # Save last sequence of the file (loop ends)
            if current_accession is not None:
                f_out.write(f"{current_accession}\t{current_length}\n")
        print(f"Success: File processed. Results saved to {output_file}.")
    
    except FileNotFoundError:
        print(f"Error: The file {input_file} was not found. Make sure that it exists in the specified path.")
                    
if __name__ == "__main__":
    # Configure arguments for command line execution to be dynamic
    # but provide default name that the lab instructions specify if no arguments are given.
    input_fasta = sys.argv[1] if len(sys.argv) > 1 else "PM_50.fasta"
    output_tsv = sys.argv[2] if len(sys.argv) > 2 else "sequences.tsv"

    parse_fasta(input_fasta, output_tsv)
                