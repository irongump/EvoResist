#1st step: preprocess the data

#2nd step: mutation calling
ls ../data/strain_ids/L*txt| while read line; do
    sbatch fastq2allres.sl $line
done
#for sra data, if you have fastq data, please run fastq2allres_cgsa.sl
#optional step
#sbtach fastq2allres_cgsa.sl

#3rd step: build lineage tree
lineages=$(ls ../data/strain_ids/L*txt| sed 's/.*\///'| sed 's/.txt//')

for lineage in $lineages; do
    echo $lineage
    sbatch bootstrap_tree.sl $lineage
done

#4th get trees' nodes' mutation information
for lineage in $lineages; do
    echo $lineage
    sbatch branch_mutation.sl $lineage
done

#5th get ancestor mutation information
sbatch ancestor_mutation.sl

#6th merge mutation information and get convergent mutation counts per genome site
sbatch merge_ann.sl

#7th run simulation to get convergent mutation null expectation distribution
mkdir -p output/simulation
cd output/simulation
sh ../scripts/simulation_GTR.sh
