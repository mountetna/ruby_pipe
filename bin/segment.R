#!/opt/R/R-latest/bin/Rscript --no-save

args = commandArgs(TRUE)

print.noquote("Segmenting...")

suppressMessages(library(DNAcopy))

# okay, using this, first generate the CBS segmentation. Return the segmentation object.
segCBS=function(cna) {
	return(segment(cna,verbose=2,alpha=0.05,nperm=10000,undo.splits='sdundo',undo.SD=1.5))
}

logr_file = args[1]
rdata_file = args[2]
seg_file = args[3]
sname = args[4]

# load coverage data (made using coverageBed from the BAMs)
tumor_logr=read.table(logr_file,header=T,as.is=T)

# clean out Y chromosomes
tumor_logr = tumor_logr[tumor_logr$chr != "chrY",]

# generate 
print.noquote("Creating CNA...")
cna=CNA(tumor_logr$logr,tumor_logr$chr,as.integer((tumor_logr$start+tumor_logr$stop)/2),"logratio",sname)
smCNA=smooth.CNA(cna)
cbs_seg=segCBS(smCNA)

save(cbs_seg,file=rdata_file)
write.table(cbs_seg$output,file=seg_file,sep="\t",eol="\n",quote=F,row.names=F)
