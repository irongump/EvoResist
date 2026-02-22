library(dplyr)
library(ggpubr)

#stat convergence for ann file
annfile <- "all_ann.txt"
ann <- read.delim(annfile,header=F,sep='\t')
names(ann) <- c('pos','ref','alt','codon_pos','type','codon_change','locus','gene','ano','grp')
ann <- ann %>% mutate(id=paste(pos,alt,sep='_'))
c_ann <- ann %>% group_by(id) %>% summarise(count=n())
c_ann <- c_ann %>% left_join(ann,multiple = 'first',by='id',) %>% 
  select(c('gene','codon_pos','type','pos','ref','alt','count'))

write.table(c_ann,'all_ann_convergent.txt',quote=F,row.names = F)
gghistogram(c_ann,x='count',binwidth = 1) + 
  scale_y_log10() + scale_x_continuous(limits=c(1,100))
