#!/usr/bin/env python3
import sys
import os

if len(sys.argv) != 3:
    sys.exit(f"usage: {sys.argv[0]} <db_mutation> <snp_dir>")

fasta_file = "../data/tb.ancestor.fasta"
with open(fasta_file, "r") as fa:
    genome = "".join(line.strip() for line in fa if not line.startswith(">"))

mutation_file = sys.argv[1]
output_snp_dir = sys.argv[2]

if not os.path.exists(output_snp_dir):
    os.makedirs(output_snp_dir)

mutations_by_sid = {}

with open(mutation_file, "r") as infile:
    for line in infile:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        sid = parts[0]
        mutation_parts = parts[2].split("_")
        pos = int(mutation_parts[0])
        alt = mutation_parts[1]
        ref = genome[pos - 1]
        mutation_line = f"{pos}\t{ref}\t{alt}\n"
        mutations_by_sid.setdefault(sid, []).append(mutation_line)

for sid, lines in mutations_by_sid.items():
    output_file = os.path.join(output_snp_dir, f"{sid}.snp")
    with open(output_file, "w") as outfile:
        outfile.writelines(lines)
