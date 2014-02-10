#!/opt/R/R-latest/bin/Rscript --no-save
source("/home/changmt/seqCBS/CAST_lib.R")
source("/home/changmt/seqCBS/sdundo.R")
sname="SeqCBS"

nor_infile=commandArgs(TRUE)[1]
tum_infile=commandArgs(TRUE)[2]
#prefix=commandArgs(TRUE)[3]

if(!file.exists(nor_infile) || !file.exists(tum_infile)) {
	logC("ERROR",sname,"Input files not found!")
	q(save="n")
}

# Load necessary items and extract the components of the file name
logC("INFO",sname,"Formatting sample and contig details")
#suppressMessages(library(seqCBS))
suppressMessages(library("clue", lib.loc="/home/changmt/seqCBS/bin"))
suppressMessages(library("seqCBS", lib.loc="/home/changmt/seqCBS/bin"))
sample_name=gsub("(^.*/|[.].*$)","",tum_infile)
#base=paste("output/",prefix,sep="")
contig=gsub("^.*[.]","",tum_infile)

# Segment
logC("INFO",sname,"Begin segmenting tumor and matched normal reads")
dn=readSeq(nor_infile,"Chiang"); # Matched normal reads
dt=readSeq(tum_infile,"Chiang"); # Tumor reads
ds=ScanCBS(dt$seqF,dn$seqF,statistic="binomial",timing=TRUE)
#ci=BayesCptCI(dt$seqF,dn$seqF,ds,verbose=TRUE)
#ScanCBSPlot(dt$seqF,dn$seqF,ds,paste(sname,contig,sep=","),"TCGA-A6-2674 chr16",ci,localSeparatePlot=FALSE,width=8,height=11.5)

# Olshen preliminary sdundo code, I adjusted the copy number value 
logC("INFO",sname,"Undoing hyper-segmentation")
all5pos=c(dt$seqF,dn$seqF)
all5label=c(rep(TRUE,length(dt$seqF)),rep(FALSE,length(dn$seqF)))
oidx=order(all5pos)
all5pos=all5pos[oidx]
all5label=all5label[oidx]
p=length(dt$seqF)/(length(dt$seqF)+length(dn$seqF))
baseline=p/(1-p)
#ds.changepoints.seqcbs=changepoints.seqcbs(ds,all5pos,all5label,0.15)
ds.changepoints.seqcbs=changepoints.seqcbs(ds,all5pos,all5label,0.15,1,baseline)

# Plot the contig and reformate output for CBS backward compatibility
logC("INFO",sname,"Plot and reformat for CBS backward compatibility")
d=NULL
d$data=list()
maxVal=max(c(dt$seqF,dn$seqF))
scale=10^6
grid.fix=seq(1,maxVal,length.out=10000)
d$data$chrom=rep(contig,length(grid.fix))
d$data$maploc=round(grid.fix)
gridSize=grid.fix[2]-grid.fix[1]
casesCountInGrid=getCountsInWindow(dt$seqF,0,maxVal,gridSize,sorted = FALSE)
controlCountInGrid=getCountsInWindow(dn$seqF,0,maxVal,gridSize,sorted = FALSE)
PInGrid=casesCountInGrid/(casesCountInGrid+controlCountInGrid)
PInGrid[is.nan(PInGrid)] = 0 
relCNInGrid=PInGrid/(1-PInGrid)/(p/(1 - p))
relCNInGrid[is.nan(relCNInGrid) | !is.finite(relCNInGrid) | relCNInGrid <= 0]=1
relCNInGrid=log2(relCNInGrid)
d$data[[sample_name]]=relCNInGrid
d$data=as.data.frame(d$data)

#pdf(file=paste(base,"__",contig,".pdf",sep=""),width=8.5,height=6)
#plot(grid.fix/scale,relCNInGrid,type="p",pch='.',cex=1,ylim=c(-2,2),ylab="log2 relative CN",xlab="Position (Mb)")
tauHat=ds$tauHat
tauHat=tauHat[c(-1,-length(tauHat))]
plotTauHatInd=c(min(c(dt$seqF,dn$seqF)),tauHat,maxVal)%/%gridSize
plotTauHatInd=sapply(plotTauHatInd,function(x) max(x,1))
plotTauHatInd=sapply(plotTauHatInd,function(x) min(x,max(grid.fix)))
plotTauHat=grid.fix[plotTauHatInd]/scale
relCN=ds$relCN
relCN[relCN<=0]=1
n=length(ds$tauHat)-1
d$output=list()
d$output$ID=rep(sample_name,n)
d$output$chrom=rep(contig,n)
d$output$loc.start=ds$tauHat[-length(ds$tauHat)]
d$output$loc.end=ds$tauHat[-1]
d$output$num.mark=rep(NA,n)
d$output$seg.mean=log2(relCN)
d$output=as.data.frame(d$output)
#lines(plotTauHat,log2(c(relCN,relCN[length(relCN)])),lwd=1.5,type="s",col="red3")

tauHat=ds.changepoints.seqcbs$taus
tauHat=tauHat[c(-1,-length(tauHat))]
plotTauHatInd=c(min(c(dt$seqF,dn$seqF)),tauHat,maxVal)%/%gridSize
plotTauHatInd=sapply(plotTauHatInd,function(x) max(x,1))
plotTauHatInd=sapply(plotTauHatInd,function(x) min(x,max(grid.fix)))
plotTauHat=grid.fix[plotTauHatInd]/scale
relCN=ds.changepoints.seqcbs$cn
relCN[relCN<=0]=1
n=length(ds.changepoints.seqcbs$taus)-1
d$undid=list()
d$undid$ID=rep(sample_name,n)
d$undid$chrom=rep(contig,n)
d$undid$loc.start=ds.changepoints.seqcbs$taus[-length(ds.changepoints.seqcbs$taus)]
d$undid$loc.end=ds.changepoints.seqcbs$taus[-1]
d$undid$num.mark=rep(NA,n)
d$undid$seg.mean=log2(relCN)
d$undid=as.data.frame(d$undid)
#lines(plotTauHat,log2(c(relCN,relCN[length(relCN)])),lwd=1.5,type="s",col="dodgerblue")
#null=dev.off()

#adding num.marks per windows defined by seqCBS in d$output
for(i in 1:nrow(d$undid)) {
  tm=d$data[which(d$data$maploc > d$undid$loc.start[i]),]
  tm=tm[which(tm$maploc < d$undid$loc.end[i]),]
  d$undid$num.mark[i]=nrow(tm)
}
#often the first row has the same start and stop position (or < 100bp
#if this the case, this is not useful and will be removed
d$undid=d$undid[-which(d$undid$num.mark==0),]

d$call="Prostate/lowpass seqCBS invokation pipeline"
class(d$data)=c("CNA","data.frame")
class(d)="DNAcopy"
#save(d,file=paste("/home/changmt/seqCBS/",base,"__",contig,".Rdata",sep=""))
save(d,file=commandArgs(TRUE)[3])
logC("INFO",sname,"Completing segmentation...")
