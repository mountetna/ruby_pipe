#!/opt/R/R-latest/bin/Rscript --no-save
suppressMessages(library(ABSOLUTE))

args = commandArgs(TRUE)

# Get input and set globals

if (length(args) < 4) {
print("Usage: absolute.R <sample_name> <seg_file> <maf_file> <results_dir>")
quit()
}

sample_name=args[1]
seg_file=args[2]
maf_file=args[3]
results_dir=args[4]

platform="Illumina"
pd="cancer"
cnt="total"

write.table(file.path(results_dir, paste(sample_name,".maf")))

maf=read.table(maf_file,header=T,as.is=T,sep="\t")
maf_fname=paste("input/",sample_name,".maf",sep="")
write.table(maf[maf$Tumor_Sample_Barcode==sample_name,],file=maf_fname,sep="\t",eol="\n",quote=F,row.names=F)

# Run ABSOLUTE 
RunAbsolute(
	seg.dat.fn=seg_file,
	maf.fn=maf_fname,
	sample.name=sample_name,
	output.fn.base=sample_name,
	results.dir=results_dir,
	primary.disease=pd,
	platform=platform,
	copy_num_type=cnt,
	sigma.p=0.01,
	max.sigma.h=0.02,
	min.mut.af=0.05,
	min.ploidy=0.95,
	max.ploidy=10,
	max.as.seg.count=15000,
	max.non.clonal=0.9,
	max.neg.genome=0.9,
	verbose=TRUE
)