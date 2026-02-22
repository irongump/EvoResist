#!/bin/bash

#SBATCH -p hov
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 3-00:00:00
#SBATCH --mem=20G

#working directory lineage_ann
cat Lineage*ann > all_ann.txt
cat 100k_5000_anc.ann >> all_ann.txt

#stat convergent times
module load r/4.4.0
Rscript ../scripts/stat_convergent.R #input all_ann.txt, output all_ann_convergent.txt

#filter low quality mutations
Rscript ../scripts/filter_low_freq_pos.R all_ann_convergent.txt all_ann_convergent_flt.txt