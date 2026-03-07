#filter out low quality samples

library(dplyr)

args <- commandArgs(T)
qfile <- args[1]
fname <- basename(qfile)
depth_file <- args[2]
depth_thres <- args[3]
fname <- gsub('.txt','',fname)
fname <- paste0(fname,'_depth_',depth_thres,'.txt')
need_remove_sample <- read.delim('data/need_remove_samples.txt',header=F)

sample <- read.delim(qfile,header=F)


depth <- read.delim(depth_file, header=T, sep='\t')
low_depth <- depth %>% filter(as.numeric(Sequencing_depth) <= as.numeric(depth_thres) | as.numeric(Genome_coverage) <= 0.9) #average depth <=30, 10x coverage <= 0.9
sample_flt <- sample %>% filter(!(V1 %in% need_remove_sample$V1)) %>% filter(!(V1 %in% low_depth$ID))

write.table(sample_flt,fname,quote=F,row.names=F,col.names = F)
