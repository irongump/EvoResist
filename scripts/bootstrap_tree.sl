#!/bin/bash
#SBATCH -p hov
#SBATCH -N 1
#SBATCH -c 64
#SBATCH -t 10-00:00:00
#SBATCH --mem=200G
#SBATCH -J iqtree_job
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# Usage: sbatch script.sh <lname> <strain_list> [<start_tree.nwk>]
# Build phylogeny using IQ-TREE, root on tb.ancestor, bootstrap=1000

set -euo pipefail   # Exit on error, report undefined variables

# -- Parameter check --
lname="${1:?Error: lname (lineage name) is required as \$1}"
strain_list="${2:?Error: strain_list is required as \$2}"

# -- Directory initialization --
output_dir="output"
tree_dir="${output_dir}/lineage_tree"
cfa_dir="${output_dir}/lineage_cfa"
snv_dir="${output_dir}/snv"
mkdir -p "$tree_dir" "$cfa_dir"

# -- Initialize temporary files (prevent appending to stale data) --
pos_file="${lname}.pos"
snp_strain_file="${lname}_snp_strain.txt"
cfa_strain_file="${lname}_cfa_strain.txt"
> "$pos_file"
> "$snp_strain_file"
> "$cfa_strain_file"

# -- Generate position file --
while IFS= read -r strain; do
    snp_file="${snv_dir}/${strain}.snp"
    if [ -f "$snp_file" ]; then
        awk '{print $1}' "$snp_file" >> "$pos_file"
        echo "$strain" >> "$snp_strain_file"
    else
        echo "Warning: ${strain} snp does not exist, skipping." >&2
    fi
done < "$strain_list"

sort -nu "$pos_file" -o "$pos_file"

# -- Generate CFA strain list --
while IFS= read -r strain; do
    cfa_file="${cfa_dir}/${strain}.cfa"
    if [ -f "$cfa_file" ]; then
        echo "$cfa_file" >> "$cfa_strain_file"
    else
        echo "Warning: ${strain} cfa does not exist, skipping." >&2
    fi
done < "$snp_strain_file"
echo "./data/tb.ancestor.concatenated.fasta" >> "$cfa_strain_file"

# -- Run snp2cfa --
module load python/3.9.6
python ./src/snp2cfa.py "$pos_file" "$cfa_strain_file"

# -- IQ-TREE phylogenetic analysis --
cfa="${lname}.fa"
if [ ! -f "$cfa" ]; then
    echo "Error: Expected CFA file $cfa not found after snp2cfa.py!" >&2
    exit 1
fi

name=$(basename "${cfa%.*}")

module load iqtree

echo "Starting IQ-TREE analysis for $name at $(date)"

iqtree2 -s "$cfa" \
        $TREE_OPT \
        -m GTR+F+R4 \
        --seqtype DNA \
        --prefix "${tree_dir}/${name}_btp" \
        --mem 190G \
        -T "${SLURM_CPUS_PER_TASK:-64}" \
        --ancestral \
        -o tb \
        -af fasta \
        -bb 1000 \
        -alrt 1000

echo "IQ-TREE finished at $(date)"