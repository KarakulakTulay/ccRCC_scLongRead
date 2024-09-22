import argparse

def parse_fasta_headers(fasta_file):
    headers = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.split()[0][1:]  # Get the first part of the header
                headers.add(header)
    return headers

def find_highest_scores(scores_file, fasta_headers):
    highest_scores = {}
    highest_score_headers = {}

    with open(scores_file, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                header = parts[0]
                if header in fasta_headers:
                    score = float(parts[3])  
                    pb_id = header.split("|")[0]

                    if pb_id not in highest_scores or score > highest_scores[pb_id]:
                        highest_scores[pb_id] = score
                        highest_score_headers[pb_id] = header

    return highest_score_headers

def extract_best_sequences(fasta_file, output_file, highest_score_headers):
    with open(fasta_file, 'r') as f, open(output_file, 'w') as out:
        write_sequence = False
        for line in f:
            if line.startswith('>'):
                header = line.split()[0][1:]  # Get the first part of the header
                pb_id = header.split("|")[0]
                if pb_id in highest_score_headers and header == highest_score_headers[pb_id]:
                    write_sequence = True
                    out.write(line)  # Write the header line
                else:
                    write_sequence = False
            elif write_sequence:
                out.write(line)  # Write the sequence lines

def main():
    parser = argparse.ArgumentParser(description="Extract highest scoring sequences from a FASTA file based on scores in a separate file.")
    parser.add_argument("fasta_file", type=str, help="Path to the FASTA file.")
    parser.add_argument("scores_file", type=str, help="Path to the scores file.")
    parser.add_argument("output_file", type=str, help="Path to the output FASTA file.")

    args = parser.parse_args()

    fasta_headers = parse_fasta_headers(args.fasta_file)
    highest_score_headers = find_highest_scores(args.scores_file, fasta_headers)
    extract_best_sequences(args.fasta_file, args.output_file, highest_score_headers)

if __name__ == "__main__":
    main()