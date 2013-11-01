#!/opt/R/R-latest/bin/Rscript --no-save

args = commandArgs(TRUE)
lib_dir = args[1]
cmd = args[2]

if (length(args) < 1) {
  print("Usage:")
  print("  deseq <lib_dir> doDeseq <counts_file> <exp_name> <ctrl_name> <fdr_cutoff> <diff_file>")
  quit()
}

doDeseq=function(counts_file,exp_name,ctrl_name,fdr,diff_file) {
	counts = read.table(counts_file,header=T,sep="\t",row.names="gene_id")
	# trim 
	counts = counts[ , grep(paste(paste(exp_name,".r",sep=""),paste(ctrl_name,".r",sep=""),sep="|"),colnames(counts)) ]
	conditions = factor( c("experiment", "control")[ as.integer(grepl(ctrl_name,colnames(counts)))+1 ] )
	library(DESeq)
	cds = newCountDataSet( counts, conditions )
	cds = estimateSizeFactors( cds )
	cds = estimateDispersions( cds )
	res = nbinomTest( cds, "control", "experiment" )
	resSig = res[ which(res$padj < fdr), ]
	write.table(resSig[ order(resSig$pval), ],file=diff_file,row.names = F,sep="\t",quote=F)
}

do.call(cmd,as.list(args[-1:-2]))
