#!/opt/R/R-latest/bin/Rscript --no-save
suppressMessages(library(ABSOLUTE))

args = commandArgs(TRUE)
lib_dir = args[1]
cmd = args[2]

# Get input and set globals

if (length(args) < 1) {
  print("Usage:")
  print("  absolute.R <lib_dir> callSample <sample_name> <seg_file> <maf_file> <results_dir>")
  print("  absolute.R <lib_dir> createReview <cohort_name> <results_dir> <list of filenames>")
  quit()
}

extractReview = function(reviewed_file, userid, mod_file, review_dir, cohort_name) {
  ExtractReviewedResults(
    reviewed.pp.calls.fn=reviewed_file,
    analyst.id=userid,
    modes.fn=mod_file,
    out.dir.base=review_dir,
    obj.name=cohort_name,
    copy_num_type="total"
    )
}

createReview = function(cohort_name, results_dir, files) {

  CreateReviewObject(
    obj.name = cohort_name,
    absolute.files = files,
    indv.results.dir = results_dir,
    copy_num_type="total",
    verbose=TRUE
  )
}

callSample = function(sample_name, seg_file, maf_file, results_dir) {
  platform="Illumina"
  pd="cancer"
  cnt="total"

  #maf=read.table(maf_file,header=T,as.is=T,sep="\t")
  #maf_fname=paste("input/",sample_name,".maf",sep="")
  #write.table(maf[maf$Tumor_Sample_Barcode==sample_name,],file=maf_fname,sep="\t",eol="\n",quote=F,row.names=F)

  # Run ABSOLUTE 
  RunAbsolute(
    seg.dat.fn=seg_file,
    maf.fn=maf_file,
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
}

do.call(cmd,as.list(args[-1:-2]))
