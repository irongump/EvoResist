#!/bin/bash

#SBATCH -p hov
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 3-00:00:00
#SBATCH --mem=20G

#get ancestor mutations
python ../scripts/getrefbase_per_node.py ../data/100k_5000_anc_db_mutation2_rmanc.txt 100k_5000_anc ../data/tb.ancestor.fasta #this script will mkdir ${name} and output ${name}/perNode.snp
find "100k_5000_anc" -maxdepth 1 -name "*snp" -print0 | parallel -0 -j 24 '
    python ../scripts/remove_low_ebr.py "../../data/RLC_lowmapK50E4_H37Rv_pos.txt" "{}" > "{= s/\.snp$// =}_rle.snp" && 
    perl ../scripts/mtbc_translate/0_MTBC_Annotation_mtbc_4411532_corrected.pl "{= s/\.snp$// =}_rle.snp" > "{= s/\.snp$// =}.ann" &&
    sed -i "/^$/d" "{= s/\.snp$// =}.ann"
'
cat 100k_5000_anc/*.ann > 100k_5000_anc.ann

#remove temp files
rm -r 100k_5000_anc