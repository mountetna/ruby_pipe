#!/usr/bin/env Rscript
sample=matrix(numeric(0), 0,0) 
args=commandArgs(trailingOnly = TRUE)
for(i in 2:length(args)) {
  load(file=args[i])
  table=d$undid
  colnames(table)=c('ID','Chromosome','Start','End','Num_Probes','Segment_Mean')
  sample=rbind(sample,table)
}
write.table(sample,file=args[1],sep="\t",row.names=F,quote=F)


