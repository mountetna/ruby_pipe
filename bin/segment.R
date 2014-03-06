#!/opt/R/R-latest/bin/Rscript --no-save

args = commandArgs(TRUE)
lib_dir = args[1]
cmd = args[2]

if (length(args) < 1) {
  print("Usage:")
  print("  segment.R <lib_dir> doSegCbs <logr_file> <rdata_file> <seg_file> <sample_name>")
  quit()
}

suppressMessages(library(DNAcopy))
suppressMessages(library(PSCBS))
source(paste(lib_dir,"bin", "chrom.R",sep="/"))
# okay, using this, first generate the CBS segmentation. Return the segmentation object.
segCBS=function(cna) {
	return(segment(cna,verbose=2,alpha=0.05,nperm=10000,undo.splits='sdundo',undo.SD=5))
}

segPSCBS=function(tumor_logr,tumor_baf,normal_baf) {
	print.noquote("Running PSCBS.")
        verbose <- Arguments$getVerbose(-10*interactive(), timestamp=TRUE)
	# get the baf and CNR into one matrix
	return(segmentByPairedPSCBS(
		CT=(3*tumor_logr$ratio),
		betaT=tumor_baf$BAF, 
		betaN=normal_baf$BAF, 
		#muN=rep(0.5,length(tumor_logr$ratio)),
		chromosome=chr_int(tumor_baf$chromosome), 
		x=tumor_baf$position, 
		#tbn=F,
		alphaTCN=0.009, alphaDH=0.001, undoTCN=25, undoDH=25,
		#..., 
		#flavor=c("tcn&dh", "tcn,dh", "sqrt(tcn),dh", "sqrt(tcn)&dh", "tcn"),
		flavor="tcn&dh",
                joinSegments=T,
                verbose=verbose
		))
}

load_logr_file = function(logr_file) {
	tumor_logr=read.table(logr_file,header=T,as.is=T,sep="\t")
	tumor_logr = tumor_logr[tumor_logr$chr != "chrY",]
	return(tumor_logr)
}

smoothCNA=function(tumor_logr, sample_name) {
	return(smooth.CNA(
		CNA(tumor_logr$logr, tumor_logr$chr, as.integer((tumor_logr$start+tumor_logr$stop)/2), "logratio", sample_name)
	))
}

# fix sorting of chromosomes
sort_cbs_seg = function(cbs_seg) {
	library(gtools)
	cbs_seg$data = cbs_seg$data[ mixedorder(cbs_seg$data$chrom), ]
	cbs_seg$output = cbs_seg$output[ mixedorder(cbs_seg$output$chrom), ]
	return(cbs_seg)
}

recenter_cbs_seg = function(cbs_seg) {
	z = density(cbs_seg$output$seg.mean,bw="SJ")
	cbs_seg$output$seg.mean = cbs_seg$output$seg.mean - z$x[z$y == max(z$y)]
	return(cbs_seg)
}

save_seg_obj = function(cbs_seg,seg_file) {
	seg_obj = cbs_seg$output
	colnames(seg_obj) = c("ID","Chromosome", "Start","End","Num_Probes","Segment_Mean")
	write.table(seg_obj,file=seg_file,sep="\t",eol="\n",quote=F,row.names=F)
}

# generate a CBS segmentation from a logr file and save it
doSegCbs=function(logr_file, rdata_file, seg_file, sample_name) {
	print.noquote("Segmenting...")

	# load coverage data (made using coverageBed from the BAMs)
	tumor_logr = load_logr_file(logr_file)
	# generate 
	print.noquote("Creating CNA...")
	cbs_seg=sort_cbs_seg(segCBS( smoothCNA(tumor_logr, sample_name) ))

	# recenter the chromosome segments by estimating the mode
	cbs_seg=recenter_cbs_seg(cbs_seg)

	save(cbs_seg,file=rdata_file)

	save_seg_obj(cbs_seg,seg_file)
}

doSegPscbs=function(logr_file,tumor_baf_file,normal_baf_file,rdata_file) {
	print.noquote("Segmenting with PSCBS...")

	tumor_logr = load_logr_file(logr_file)

	tumor_baf = load_baf(tumor_baf_file)
	normal_baf = load_baf(normal_baf_file)

	pscbs_seg = segPSCBS(tumor_logr, tumor_baf, normal_baf)

    save(pscbs_seg,file=rdata_file)
#        deltaAB=estimateDeltaAB(pscbs_seg,flavor="qq(DH)")
#        fit=callAB(pscbs_seg,delta=deltaAB,verbose=verbose)
#        save(fit,file=paste(rdata_file,"_fit",sep=""))

#       deltaLOH=estimateDeltaLOH(pscbs_seg, flavor="minC1|nonAB")
#        fit=callLOH(pscbs_seg,delta=deltaLOH,verbose=verbose)
#        save(fit,file=paste(rdata_file,"LOH",sep="_"))

#        tab=getSegments(pscbs_seg,simplify=T)
#        write.table(tab,file=paste(rdata_file,".txt",sep=""),row.names=F,quote=F)

#        chrTag <- sprintf("Chr%s", seqToHumanReadable(getChromosomes(pscbs_seg)));
#        toPNG(paste(rdata_file,"png",sep="."), tags=c(chrTag, "PairedPSCBS"), width=840, aspectRatio=0.6, {
#           plotTracks(pscbs_seg);
#        });
}

getASCAT=function(ascatpcf) {
	ascatoutput = ascat.runAscat(ascatpcf)
	return(ascatoutput)
}

ascat_dir=function() {

}

segASCAT=function(tumor_ba_logr,tumor_baf,normal_ba_logr,normal_baf) {
	source("/taylorlab/scripts/ruby_pipe/bin/ASCAT2.1/ascat.R")
	bc = ascat.loadData(tumor_ba_logr,tumor_baf,normal_ba_logr,normal_baf)
	pcf = ascat.aspcf(bc)
	return(pcf)
}

make_ba_logr=function(logr,baf) {
	baf_to_logr = map_to_logr(baf,logr)
	ba_logr = list(chromosome=gsub("chr","",baf$chromosome),
		  position=baf$position,
		  logR = NA)
	ba_logr$logR[1:length(baf$position)] = NA
	ba_logr$logR[baf_to_logr[,1]] = logr$logr[baf_to_logr[,2]]
	ba_logr = as.data.frame(ba_logr)
	return(ba_logr)
}

doAscatPurityPloidy=function(tumor_baf_file,normal_baf_file,tumor_logr_file,normal_logr_file,exon_intervals,
		ascat_rdata_file,ascat_txt) {

	tumor_baf = load_baf(tumor_baf_file)
	normal_baf = load_baf(normal_baf_file)
	tumor_logr = load_logr_file(tumor_logr_file)
	normal_logr = load_logr_file(normal_logr_file)
	
	normal_ba_logr = make_ba_logr(normal_logr, normal_baf)
	tumor_ba_logr = make_ba_logr(tumor_logr, tumor_baf)

	ascat_seg = segASCAT(tumor_ba_logr,tumor_baf,normal_ba_logr,normal_baf)
	ascat_out = getASCAT(ascat_seg)
	save(ascat_out,file=ascat_rdata_file)
	write.table(list(tumor_frac=ascat_out$aberrantcellfraction,
			ploidy=ascat_out$ploidy),
			file=ascat_txt,
			sep="\t",quote=F,row.names=F)
}

load_baf = function(baf_file) {
	return(read.table(baf_file,header=T,as.is=T,sep="\t"))
}

map_to_exons = function(baf,exon_intervals) {
	suppressMessages(library(GenomicRanges))

	exons = read.table(exon_intervals,header=F,as.is=T,sep="\t")[ , 1:3]
	colnames(exons)=c("chrom","start","end")
	return(as.matrix(
		findOverlaps(
			# tumor range
			GRanges(seqnames=baf[,1],ranges=IRanges(baf[,2],baf[,2])),
			# exon range
			GRanges(seqnames=exons[,1],ranges=IRanges(exons[,2],exons[,3]))
		)
	))
}

map_to_logr = function(baf,logr) {
	suppressMessages(library(GenomicRanges))

	return(as.matrix(
		findOverlaps(
			# tumor range
			GRanges(seqnames=baf[,"chromosome"],ranges=IRanges(baf[,"position"],baf[,"position"])),
			# exon range
			GRanges(seqnames=logr[,"chr"],ranges=IRanges(logr[,"start"],logr[,"stop"]))
		)
	))
}

add_gc = function(cov, gc) {
  cov$pctGC=gc$X5_pct_gc
  cutoff=mean(cov$reads)+1.5*sd(cov$reads)
  cov=cov[-which(cov$reads>cutoff),]
  cov=cov[-which(cov$pctGC < 0.33),]
  return(cov)
}

gc_correct=function(coverage,span=0.75,degree=2) {
  lo=loess(coverage$reads ~ coverage$pctGC,
  	span=span,
	degree=degree,
	data.frame(reads=coverage$reads, gc=coverage$pctGC))
  center=mean(coverage$reads)
  adj=lo$fitted-center
  coverage$reads=coverage$reads-adj
  return(coverage)
}

doGcCorrect=function(cov_file, corr_file, gc_file) {
	cov = read.table(cov_file,header=F,sep="\t")
	colnames(cov) = c("contig","start","end","reads")
	gc = read.table(gc_file,header=T,sep="\t")
	cov = add_gc(cov, gc)
	cov = gc_correct(cov)
	write.table(cov, corr_file, quote=F, sep="\t",col.names = FALSE, row.names = FALSE)
}

do.call(cmd,as.list(args[-1:-2]))
