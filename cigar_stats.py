import sys
import re

def main():
    # Setup command line arguments with defaults for Snakemake compatibility
    input_alignment = sys.argv[1] if len(sys.argv) > 1 else "alignments.tsv"
    out_stats = sys.argv[2] if len(sys.argv) > 2 else "cigar_stats.tsv"
    out_report = sys.argv[3] if len(sys.argv) > 3 else "report.txt"

    # Dictionaries to track the sequences for our final report
    min_dist_record = None
    max_dist_record = None
    max_match_record = None

    try:
        with open(input_alignment, 'r') as f_in, open(out_stats, 'w') as f_stats:
            # Write header for the stats file
            f_stats.write("Accession\tM\tI\tD\tX\n")

            # Skip the header line of alignments.tsv
            next(f_in)

            for line in f_in:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue

                accession = parts[0]
                distance = int(parts[1])
                cigar = parts[2]

                # Extract counts of M, I, D, X from the CIGAR string with regex
                m_count = len(re.findall(r'M', cigar))
                i_count = len(re.findall(r'I', cigar))
                d_count = len(re.findall(r'D', cigar))
                x_count = len(re.findall(r'X', cigar))

                # Write to cigar_stats.tsv
                f_stats.write(f"{accession}\t{m_count}\t{i_count}\t{d_count}\t{x_count}\n")

                # Store the record for the report
                record = {
                    'accession': accession,
                    'distance': distance,
                    'cigar': cigar,
                    'm': m_count,
                    'i': i_count,
                    'd': d_count,
                    'x': x_count
                }

                # Update our trackers for the smallest/largest metrics
                if min_dist_record is None or distance < min_dist_record['distance']:
                    min_dist_record = record
                if max_dist_record is None or distance > max_dist_record['distance']:
                    max_dist_record = record
                if max_match_record is None or m_count > max_match_record['m']:
                    max_match_record = record

        # Generate the summary report.txt
        with open(out_report, 'w') as f_rep:
            def write_record(title, rec):
                f_rep.write(f"{title}:\n")
                f_rep.write(f"{rec['accession']}\n")
                f_rep.write(f"Distance: {rec['distance']}\n")
                f_rep.write(f"Matches: {rec['m']}\n")
                f_rep.write(f"Insertions: {rec['i']}\n")
                f_rep.write(f"Deletions: {rec['d']}\n")
                f_rep.write(f"Substitutions: {rec['x']}\n\n")

            write_record("Sequence with smallest edit distance", min_dist_record)
            write_record("Sequence with largest edit distance", max_dist_record)
            write_record("Sequence with largest number of matches", max_match_record)

        print(f"Success: Regular Expression analysis completed.")
        print(f" - Stats saved to: {out_stats}")
        print(f" - Report saved to: {out_report}")

    except FileNotFoundError:
        print(f"Error: The file '{input_alignment}' was not found.")

if __name__ == "__main__":
    main()