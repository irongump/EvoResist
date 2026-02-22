
#optional step: get mutation rate using 100k strain data, excluding drug resistance mutation
#get fourfold degenearate site mutation rate using 100K snp ano excluding drug resistance mutation
#python ../scripts/fourfold_dgr_rate.py 

#simulate convergent mutation counts per genome site using GTR + gamma model plus fourfold degenearate site mutation rate
python ../scripts/simulation_GTR_gamma.py 

#multiple test adjust
conda activate r4
Rscript ../scripts/fdr4simulation.R
#plot expected null convergent mutation counts bar plot
Rscript ../scripts/plot_GTR_simulation.R
