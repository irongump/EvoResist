import sys
from Bio import SeqIO

def filter_and_linearize_fasta(strain_list_path, input_fasta_path, output_fasta_path):
    """
    Reads a list of IDs, filters an input FASTA for matching records,
    and writes them to a new file where each sequence occupies a single line.
    """
    # 1. Load target strain IDs into a set for O(1) lookup speed
    target_strains = set()
    try:
        with open(strain_list_path, 'r') as f:
            for line in f:
                identifier = line.strip()
                if identifier:
                    target_strains.add(identifier)
    except FileNotFoundError:
        print(f"Error: Strain list file '{strain_list_path}' not found.")
        return

    # 2. Parse input and write matching records directly to the output
    count = 0
    try:
        with open(output_fasta_path, 'w') as output_handle:
            # SeqIO.parse reads the file record-by-record (memory efficient)
            for record in SeqIO.parse(input_fasta_path, "fasta"):
                if record.id in target_strains:
                    # 'fasta-2line' forces the sequence to stay on one line
                    SeqIO.write(record, output_handle, "fasta-2line")
                    count += 1
        
        print(f"Done. Successfully wrote {count} sequences to {output_fasta_path}")

    except Exception as e:
        print(f"An error occurred during processing: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <strain_list.txt> <input.fasta> <output.fasta>")
    else:
        filter_and_linearize_fasta(sys.argv[1], sys.argv[2], sys.argv[3])