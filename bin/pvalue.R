#!/opt/R/R-latest/bin/Rscript --no-save

library(edgeR)

cohort_file=function(cohort,suffix) {
	return(paste(cohort,suffix,sep="."))
}
output_file=function(prefix,cohort,suffix) {
	return(paste(prefix,"output",cohort_file(cohort,suffix),sep="/"))
}
metrics_file=function(prefix,cohort,suffix) {
	return(paste(prefix,"metrics",cohort_file(cohort,suffix),sep="/"))
}

load_reads=function(read_file) {
	table=read.table(read_file,header=T,as.is=T,sep="\t")
	table=table[,c(1:3,grep("rna",colnames(table)))]
	return(table)
}
pvalue_keepers=function(prefix, cohort, keepers_file) {
	if(is.na(prefix)) {
	  cat("Must specify a unique prefix string!\n")
	  q(save="no")
	}

	pseudo=0.01

	metrics=read.table(metrics_file(prefix,cohort,"qc_summary"),header=T,as.is=T,row.names=1,sep="\t")
	genes=load_reads(output_file(prefix,cohort,"normal_cov"))
	rands=load_reads(output_file(prefix,cohort,"null_cov"))

	pct_spl=as.numeric(metrics["splice_counts",])/as.numeric(metrics["total",])
	names(pct_spl)=colnames(metrics)

	idx=4:ncol(genes)
	rem=which(rands[,idx]==0)
	tmm=calcNormFactors(data.matrix(genes[,idx]),refColumn=1)
	csum=colSums(genes[,idx])
	eff_lib_size=csum*tmm
	rlibsize=(csum-(csum*pct_spl[match(names(csum),names(pct_spl))]))*tmm

	grpkm=do.call(cbind,lapply(idx,function(x) ((1e9)*genes[,x])/(eff_lib_size[x-3]*as.numeric(genes$size))))
	rrpkm=do.call(cbind,lapply(idx,function(x) ((1e9)*(rands[,x]+pseudo))/((rlibsize[x-3]+pseudo)*as.numeric(rands$size))))
	colnames(grpkm)=colnames(rrpkm)=colnames(genes)[idx]
	rownames(grpkm)=genes$gene
	rownames(rrpkm)=rands$gene
	gg=log10(as.vector(grpkm))
	rr=log10(as.vector(rrpkm))

	rcdf=ecdf(rr[-rem])
	pval1=apply(grpkm,2,function(ge) 1-rcdf(log10(ge)))
	colnames(pval1)=colnames(grpkm)
	rownames(pval1)=rownames(grpkm)

	write.table(rownames(pval1)[rowSums(pval1<0.05)==ncol(pval1)],file=keepers_file,sep="\n",quote=F,col.names=F,row.names=F)
}
