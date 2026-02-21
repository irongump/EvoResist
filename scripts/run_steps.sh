#1st step: preprocess the data

#2nd step: mutation calling
sbatch fastq2allres.sl

#3rd step: build lineage tree
sbatch bootstrap_tree.sl
