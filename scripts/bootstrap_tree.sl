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

set -euo pipefail   # 遇到错误立即退出，未定义变量报错

# ── 参数检查 ──────────────────────────────────────────────────
lname="${1:?Error: lname (lineage name) is required as \$1}"
strain_list="${2:?Error: strain_list is required as \$2}"

# ── 目录初始化 ────────────────────────────────────────────────
output_dir="output"
tree_dir="${output_dir}/lineage_tree"
cfa_dir="${output_dir}/lineage_cfa"
mkdir -p "$tree_dir" "$cfa_dir"   # mkdir -p 不需要先判断

# ── 临时文件初始化（防止追加污染）─────────────────────────────
pos_file="${lname}.pos"
snp_strain_file="${lname}_snp_strain.txt"
cfa_strain_file="${lname}_cfa_strain.txt"
> "$pos_file"           # 清空或新建
> "$snp_strain_file"
> "$cfa_strain_file"

# ── 生成 position 文件 ────────────────────────────────────────
while IFS= read -r strain; do
    snp_file="./snv/${strain}.snp"
    if [ -f "$snp_file" ]; then
        awk '{print $1}' "$snp_file" >> "$pos_file"
        echo "$strain" >> "$snp_strain_file"
    else
        echo "Warning: ${strain} snp does not exist, skipping." >&2
    fi
done < "$strain_list"   # 避免子 Shell

sort -nu "$pos_file" -o "$pos_file"   # -u 直接去重，不需要 uniq

# ── 生成 cfa strain 列表 ──────────────────────────────────────
while IFS= read -r strain; do
    cfa_file="${cfa_dir}/${strain}.cfa"
    if [ -f "$cfa_file" ]; then
        echo "$cfa_file" >> "$cfa_strain_file"
    else
        echo "Warning: ${strain} cfa does not exist, skipping." >&2
    fi
done < "$snp_strain_file"
echo "./data/tb.ancestor.concatenated.fasta" >> "$cfa_strain_file"

# ── 运行 snp2cfa ──────────────────────────────────────────────
module load python/3.9.6
python ./src/snp2cfa.py "$pos_file" "$cfa_strain_file"

# ── IQ-TREE 建树 ──────────────────────────────────────────────
cfa="${lname}.fa"   # snp2cfa.py 输出固定为 ${lname}.fa
if [ ! -f "$cfa" ]; then
    echo "Error: Expected CFA file $cfa not found after snp2cfa.py!" >&2
    exit 1
fi

name=$(basename "${cfa%.*}")

module load iqtree

echo "Starting IQ-TREE analysis for $name at $(date)"

iqtree2 -s "$cfa" \
        $TREE_OPT \
        -m GTR+F+R4 \       # 显式指定 rate categories
        --seqtype DNA \
        --prefix "${tree_dir}/${name}_btp" \   # 输出到 tree_dir
        --mem 190G \
        -T "${SLURM_CPUS_PER_TASK:-64}" \      # 防止变量未定义
        --ancestral \
        -o tb \
        -af fasta \
        -bb 1000 \
        -alrt 1000

echo "IQ-TREE finished at $(date)"