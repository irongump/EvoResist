#!/bin/bash

#SBATCH -p hov
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 3-00:00:00
#SBATCH --mem=20G

name=$1 #lineage name

#working dir output/lineage_tree
#extract ancestral seq for each node
module purge
module load anaconda
conda activate bioperl
perl ../scripts/nodes_base_locus_iqtree.pl ${name}.treefile ../lineage_cfa/${name}_delete.pos ${name}.state $cfa ${name}.db ${name}_homoplasy.txt

#extract mutations in the db file
perl ../scripts/db2mutation.pl ${name}.db > ${name}_db_mutation.txt
#paste number of leafs for each node to db_mutation file
python ../scripts/node_leafs.py ${name}.treefile ${name}_db_mutation.txt > ${name}_db_mutation2.txt

#remove ancestral mutations and get lineage specific mutations
python ../scripts/filter_lineage_defining.py ${name}_db_mutation2.txt ${name} > ${name}_db_mutation2_rmanc.txt

python ../../src/getrefbase_per_node.py ${name}_db_mutation2_rmanc.txt ${name} #this script will mkdir ${name} and output ${name}/perNode.snp
find "${name}" -maxdepth 1 -name "*snp" -print0 | parallel -0 -j 24 '
    python /proj/qliulab/MTB_phy_db/src/remove_low_ebr.py "{}" > "{= s/\.snp$// =}_rle.snp" && 
    perl /proj/qliulab/MTB_phy_db/src/mtbc_translate/0_MTBC_Annotation_mtbc_4411532_corrected.pl "{= s/\.snp$// =}_rle.snp" > "{= s/\.snp$// =}.ann" &&
    sed -i "/^$/d" "{= s/\.snp$// =}.ann"
'
cat "${name}"/*.ann > "${name}.ann"


# perl ../scripts/getrefbase.pl ${name}_db_mutation2_rmanc.txt > ${name}.snp
# python ../scripts/remove_low_ebr.py ${name}.snp > ${name}_rle.snp 
# perl ../scripts/mtbc_translate/0_MTBC_Annotation_mtbc_4411532_corrected.pl ${name}_rle.snp > {}_rmanc.ann 

#remove temp files
rm ${name}_db_mutation.txt ${name}_db_mutation2.txt
rm -r ${name}
