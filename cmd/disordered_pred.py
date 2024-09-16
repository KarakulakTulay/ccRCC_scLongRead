import os
import subprocess
import argparse

def split_and_process(file_path, output_dir, iupred_path):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(file_path, 'r') as file:
        seq = ''
        header = ''
        for line in file:
            if line.startswith('>'):
                if seq:
                    process_sequence(header, seq, iupred_path, output_dir)
                    seq = ''
                header = line.strip()
            else:
                seq += line.strip()
        if seq:
            process_sequence(header, seq, iupred_path, output_dir)



def process_sequence(header, sequence, iupred_path, output_dir):
        # Extracting the PB_id from the header
        
        pb_id = header.split('|')[0].strip().lstrip('>')
                
                    # Constructing the file name
        seq_file_name = f"{pb_id}_iupred.pep"
        seq_file_path = os.path.join(output_dir, seq_file_name)

        with open(seq_file_path, 'w') as seq_file:
            
            seq_file.write(f"{header}\n{sequence}")                                    
                
        run_iupred(seq_file_path, iupred_path, output_dir)



        os.remove(seq_file_path)



def run_iupred(seq_file_path, iupred_path, output_dir):
        # Get the base name of the sequence file (without directory path)
        
        base_name = os.path.basename(seq_file_path)

         # Replace the .pep extension with _output.txt to create the output file name
        output_file_name = base_name.replace('.pep', '_output.txt')

         # Construct the full path for the output file
        output_file_path = os.path.join(output_dir, output_file_name)

        command = f"python3 {iupred_path} -a {seq_file_path} long > {output_file_path}"
        subprocess.run(command, shell=True)

        # Print statement to confirm that the command was executed
        print(f"iupred2a run on {seq_file_path}, output saved to {output_file_path}")



def main():
    parser = argparse.ArgumentParser(description="Process pep files with iupred2a.")
    parser.add_argument("pep_file_path", help="Path to the pep file.")
    parser.add_argument("iupred_path", help="Path to the iupred2a.py script.")
    parser.add_argument("output_dir", help="Directory to store output files.")
    args = parser.parse_args()

    split_and_process(args.pep_file_path, args.output_dir, args.iupred_path)

if __name__ == "__main__":
    main()