#!/usr/bin/env Rscript

args = commandArgs(TRUE)
lib_dir = args[1]
cmd = args[2]

if (length(args) < 1) {
  print("Usage:")
  print("  deseq <lib_dir> doDeseq <counts_file> <exp_name> <ctrl_name> <fdr_cutoff> <diff_file>")
  quit()
}

replicate_name = function(name) {
	return(paste( name,"\\.r", sep=""))
}

doDeseq=function(counts_file,exp_name,ctrl_name,fdr,diff_file) {
	counts = read.table(counts_file,header=T,sep="\t",row.names="gene_id")
	# trim 
	counts = counts[ , grep(paste(replicate_name(exp_name),replicate_name(ctrl_name),sep="|"),colnames(counts)) ]
	conditions = factor( c("experiment", "control", "other")[ as.integer(grepl(ctrl_name,colnames(counts)))+1 ] )
	library(DESeq)
	cds = newCountDataSet( counts, conditions )
	cds = estimateSizeFactors( cds )
	cds = estimateDispersions( cds )
	res = nbinomTest( cds, "control", "experiment" )
	fdr = as.numeric(fdr)
	resSig = res[ which(res$padj < fdr), ]
	write.table(resSig[ order(resSig$pval), ],file=diff_file,row.names = F,sep="\t",quote=F)
}

do.call(cmd,as.list(args[-1:-2]))
