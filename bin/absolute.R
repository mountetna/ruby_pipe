#!/opt/R/R-latest/bin/Rscript --no-save
suppressMessages(library(ABSOLUTE))

args = commandArgs(TRUE)
lib_dir = args[1]
cmd = args[2]

# Get input and set globals

if (length(args) < 1) {
  print("Usage:")
  print("  absolute.R <lib_dir> callSample <sample_name> <seg_file> <results_dir>")
  print("  absolute.R <lib_dir> createReview <cohort_name> <results_dir> <list of filenames>")
  quit()
}

extractReview = function(args) {
  if (length(args) < 5) {
    print("Usage:")
    print("  absolute.R extractReview <reviewed_file> <userid> <mode_file> <review_dir> <cohort_name>")
    quit()
  }
  reviewed_file = args[1]
  userid = args[2]
  mod_file = args[3]
  review_dir = args[4]
  cohort_name = args[5]

  ExtractReviewedResults(
    reviewed.pp.calls.fn=reviewed_file,
    analyst.id=userid,
    modes.fn=mod_file,
    out.dir.base=review_dir,
    obj.name=cohort_name,
    copy_num_type="total"
    )
}

createReview = function(args) {
  if (length(args) < 3) {
  print("Usage: absolute.R createReview <cohort_name> <results_dir> [<list of files>]")
  quit()
  }

  cohort_name=args[1]
  results_dir=args[2]
  files=args[-(1:2)]

  CreateReviewObject(
    obj.name = cohort_name,
    absolute.files = files,
    indv.results.dir = results_dir,
    copy_num_type="total",
    verbose=TRUE
  )
}

callSample = function(args) {
  if (length(args) < 4) {
  print("Usage: absolute.R callSample <sample_name> <seg_file> <maf_file> <results_dir>")
  quit()
  }

  sample_name=args[1]
  seg_file=args[2]
  maf_file=args[3]
  results_dir=args[4]
  #maf_file=args[4]

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

do.call(cmd,list(args[-1,-2]))
