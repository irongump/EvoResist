import sys
from Bio import SeqIO

# Ensure correct usage
if len(sys.argv) < 4:
    print("Usage: python script.py <strain_list> <input_fasta> <output_fasta>")
    sys.exit(1)

strain_file = sys.argv[1]
old_cfa = sys.argv[2]
output_cfa = sys.argv[3]

# 1. Read the list of target strains into a set for O(1) lookup speed
strains_to_keep = set()
with open(strain_file, 'r') as f:
    for line in f:
        identifier = line.strip()
        if identifier:
            strains_to_keep.add(identifier)

# 2. Parse the old_cfa and filter sequences
# We use a generator expression to save memory for large genomic files
selected_records = []

# SeqIO.parse returns an iterator of SeqRecord objects
for record in SeqIO.parse(old_cfa, "fasta"):
    # Check if the record ID or name matches your list
    if record.id in strains_to_keep:
        selected_records.append(record)
        # Optional: remove from set to track if any strains were missing
        # strains_to_keep.remove(record.id)

# 3. Output the selected records to the new file
if selected_records:
    with open(output_cfa, "w") as output_handle:
        SeqIO.write(selected_records, output_handle, "fasta")
    print(f"Successfully wrote {len(selected_records)} sequences to {output_cfa}")
else:
    print("No matching strains found in the input FASTA file.")