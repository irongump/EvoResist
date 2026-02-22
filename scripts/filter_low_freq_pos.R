library(dplyr)
library(tidyr)
library(stringr)

args <- commandArgs(T)
input_file <- args[1] #e.g. all_ann_convergent.txt
output_file <- args[2] 

pos_freqs <- read.delim("../data/convg_5_mutation_snp_freq.txt",header=T,sep=' ') #this file is from ../convergent_check
pos_freqs$mean_freq <- as.numeric(pos_freqs$mean_freq)
pos_freqs$depth_ratio <- as.numeric(pos_freqs$depth_ratio)

#filter mutations with SNP freq < 90% or mean depth < 16
freq_90 <- pos_freqs %>% filter(mean_freq<90 | depth_ratio < 0.5 |depth_ratio > 2)
#read repeat and mobile element regions
repeat_region <- read.delim('repeat_region.txt',header=F,sep=' ')
repeat_region <- repeat_region %>% select(V7) %>% separate(V7,into=c('start','end'))
repeat_region <- repeat_region %>% mutate(start=as.numeric(start) - 10, end = as.numeric(end) + 10)
mobile_element <- read.delim('mobile_element.txt',header=F,sep=' ')
mobile_element <- mobile_element %>% select(V8) %>% separate(V8,into=c('start','end'))
mobile_element <- mobile_element %>% mutate(start=as.numeric(start) - 10, end = as.numeric(end) + 10)

repeat_mobile <- rbind(repeat_region,mobile_element)
#filter low freq
raw_ann <- read.table(input_file,header=T,sep="")
flt_ann <- raw_ann %>% filter(!(pos %in% freq_90$position)) %>% filter(gene!='Rv0007')

#filter position in repeat or mobile
library(purrr)
is_in_mobile_region <- function(pos, mobile_df) {
  any(pos >= mobile_df$start & pos <= mobile_df$end)
}

flt_ann2 <- flt_ann %>% mutate(pos=as.numeric(pos)) %>% 
  filter(!map_lgl(pos, ~is_in_mobile_region(.x, repeat_mobile)))

flt_ann2 <- flt_ann2 %>% filter(!(gene=='Rv2652c-Rv2653c')) %>% filter(pos!=1480945)
flt_ann3 <- flt_ann %>% filter(str_detect(pos,'-'))
flt_out <- rbind(flt_ann2,flt_ann3)
flt_out <- flt_out %>% select(-ancestor_node,-L5, -L6, -L7)
write.table(flt_out,output_file,quote=F,row.names = F,sep='\t')


