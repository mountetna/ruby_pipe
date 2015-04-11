
stratton_plot = function(maf_file) {
	##load in maf file
	df <- read.table(file=maf_file, header=TRUE, sep="\t",comment.char="#")
	##load in exome and whole genome context lookup files
	ex_lookup <- read.csv(file="/taylorlab/home/chapmanjs/tricount/exon_triCount.txt", header=TRUE, sep="\t")

	##Load in libraries
	library(sqldf)
	library(ggplot2)
	library(plyr)
	library(scales)

	##MUNGE DATASET

	##create dataframe with desired variables and remove NA contexts
	df = df[ which(!is.na(df$context) & nchar(as.character(df$context)) > 0), ]
	df_ctxt <- as.data.frame(table(df$context),stringsAsFactors=FALSE)
	colnames(df_ctxt) <- c("context", "Freq")
	df_ctxt = df_ctxt[ 2:nrow(df_ctxt), ]

	mods = c("A[C/A]A", "A[C/A]C", "A[C/A]G", "A[C/A]T", "C[C/A]A", "C[C/A]C", "C[C/A]G", "C[C/A]T", "G[C/A]A",
	"G[C/A]C", "G[C/A]G", "G[C/A]T", "T[C/A]A", "T[C/A]C", "T[C/A]G", "T[C/A]T", "A[C/G]A", "A[C/G]C", "A[C/G]G", "A[C/G]T", "C[C/G]A", "C[C/G]C", "C[C/G]G", "C[C/G]T",
	"G[C/G]A", "G[C/G]C", "G[C/G]G", "G[C/G]T", "T[C/G]A", "T[C/G]C", "T[C/G]G", "T[C/G]T", "A[C/T]A", "A[C/T]C", "A[C/T]G", "A[C/T]T", "C[C/T]A", "C[C/T]C", "C[C/T]G",
	"C[C/T]T", "G[C/T]A", "G[C/T]C", "G[C/T]G", "G[C/T]T", "T[C/T]A", "T[C/T]C", "T[C/T]G", "T[C/T]T", "A[T/A]A", "A[T/A]C", "A[T/A]G", "A[T/A]T", "C[T/A]A", "C[T/A]C",
	"C[T/A]G", "C[T/A]T", "G[T/A]A", "G[T/A]C", "G[T/A]G", "G[T/A]T", "T[T/A]A", "T[T/A]C", "T[T/A]G", "T[T/A]T", "A[T/C]A", "A[T/C]C", "A[T/C]G", "A[T/C]T", "C[T/C]A",
	"C[T/C]C", "C[T/C]G", "C[T/C]T", "G[T/C]A", "G[T/C]C", "G[T/C]G", "G[T/C]T", "T[T/C]A", "T[T/C]C", "T[T/C]G", "T[T/C]T", "A[T/G]A", "A[T/G]C", "A[T/G]G", "A[T/G]T",
	"C[T/G]A", "C[T/G]C", "C[T/G]G", "C[T/G]T", "G[T/G]A", "G[T/G]C", "G[T/G]G", "G[T/G]T", "T[T/G]A", "T[T/G]C", "T[T/G]G", "T[T/G]T")

	for (i in 1:length(mods)) {
		if (!mods[i] %in% df_ctxt$context) {
			new_row=c(mods[i],0)
			df_ctxt = rbind(df_ctxt,new_row)
		}
	}
	##Remove "NA" contexts - use table(df_ctxt$context) to figure out which NA categories must be removed.
	#df_ctxt <- df_ctxt[ which(df_ctxt$context
	#df_ctxt <- subset(df_ctxt, !context %in% c("N[C/T]N", "NA[C/A]NA", "NA[C/G]NA", "NA[C/T]NA", "NA[T/C]NA"))

	##Create new variables
	df_ctxt$Context_mod <- gsub("[[:punct:]]", "", df_ctxt$context)
	df_ctxt$Context_mod <- paste(substr(df_ctxt$Context_mod, 1, 2), substr(df_ctxt$Context_mod, 4, 5), sep='')


	##Append WG and exome lookup counts to df_ctxt
	df_ctxt$ctxt_count <- ex_lookup$count[match(df_ctxt$Context_mod, ex_lookup$tri)]

	##Create new variable, mutct, based on number of mutations per tumor

	##Create rescale variable which is mutrate = Freq/ctxt_count and mutrate/sum of all mutrates for individual tumors
	df_ctxt$mutrate <- as.numeric(df_ctxt$Freq)/df_ctxt$ctxt_count
	df_ctxt$norm <- df_ctxt$mutrate/sum(df_ctxt$mutrate)

	##ORDERING
	df_ctxt$context <- factor(df_ctxt$context, levels=c("A[C/A]A", "A[C/A]C", "A[C/A]G", "A[C/A]T", "C[C/A]A", "C[C/A]C", "C[C/A]G", "C[C/A]T", "G[C/A]A",
	"G[C/A]C", "G[C/A]G", "G[C/A]T", "T[C/A]A", "T[C/A]C", "T[C/A]G", "T[C/A]T", "A[C/G]A", "A[C/G]C", "A[C/G]G", "A[C/G]T", "C[C/G]A", "C[C/G]C", "C[C/G]G", "C[C/G]T",
	"G[C/G]A", "G[C/G]C", "G[C/G]G", "G[C/G]T", "T[C/G]A", "T[C/G]C", "T[C/G]G", "T[C/G]T", "A[C/T]A", "A[C/T]C", "A[C/T]G", "A[C/T]T", "C[C/T]A", "C[C/T]C", "C[C/T]G",
	"C[C/T]T", "G[C/T]A", "G[C/T]C", "G[C/T]G", "G[C/T]T", "T[C/T]A", "T[C/T]C", "T[C/T]G", "T[C/T]T", "A[T/A]A", "A[T/A]C", "A[T/A]G", "A[T/A]T", "C[T/A]A", "C[T/A]C",
	"C[T/A]G", "C[T/A]T", "G[T/A]A", "G[T/A]C", "G[T/A]G", "G[T/A]T", "T[T/A]A", "T[T/A]C", "T[T/A]G", "T[T/A]T", "A[T/C]A", "A[T/C]C", "A[T/C]G", "A[T/C]T", "C[T/C]A",
	"C[T/C]C", "C[T/C]G", "C[T/C]T", "G[T/C]A", "G[T/C]C", "G[T/C]G", "G[T/C]T", "T[T/C]A", "T[T/C]C", "T[T/C]G", "T[T/C]T", "A[T/G]A", "A[T/G]C", "A[T/G]G", "A[T/G]T",
	"C[T/G]A", "C[T/G]C", "C[T/G]G", "C[T/G]T", "G[T/G]A", "G[T/G]C", "G[T/G]G", "G[T/G]T", "T[T/G]A", "T[T/G]C", "T[T/G]G", "T[T/G]T"), ordered=TRUE)


	##Revalue context to a new category for coloring a la Alexandrov Signatures
	df_ctxt$color_ctxt <- revalue(df_ctxt$context, c("A[C/A]A"="C>A", "A[C/A]C"="C>A", "A[C/A]G"="C>A", "A[C/A]T"="C>A", "C[C/A]A"="C>A", "C[C/A]C"="C>A", "C[C/A]G"="C>A", "C[C/A]T"="C>A",
	"G[C/A]A"="C>A", "G[C/A]C"="C>A", "G[C/A]G"="C>A", "G[C/A]T"="C>A", "T[C/A]A"="C>A", "T[C/A]C"="C>A", "T[C/A]G"="C>A", "T[C/A]T"="C>A", "A[C/G]A"="C>G", "A[C/G]C"="C>G", "A[C/G]G"="C>G",
	"A[C/G]T"="C>G", "C[C/G]A"="C>G", "C[C/G]C"="C>G", "C[C/G]G"="C>G", "C[C/G]T"="C>G", "G[C/G]A"="C>G", "G[C/G]C"="C>G", "G[C/G]G"="C>G", "G[C/G]T"="C>G", "T[C/G]A"="C>G", "T[C/G]C"="C>G",
	"T[C/G]G"="C>G", "T[C/G]T"="C>G", "A[C/T]A"="C>T", "A[C/T]C"="C>T", "A[C/T]G"="C>T", "A[C/T]T"="C>T", "C[C/T]A"="C>T", "C[C/T]C"="C>T", "C[C/T]G"="C>T", "C[C/T]T"="C>T", "G[C/T]A"="C>T",
	"G[C/T]C"="C>T", "G[C/T]G"="C>T", "G[C/T]T"="C>T", "T[C/T]A"="C>T", "T[C/T]C"="C>T", "T[C/T]G"="C>T", "T[C/T]T"="C>T", "A[T/A]A"="T>A", "A[T/A]C"="T>A", "A[T/A]G"="T>A", "A[T/A]T"="T>A",
	"C[T/A]A"="T>A", "C[T/A]C"="T>A", "C[T/A]G"="T>A", "C[T/A]T"="T>A", "G[T/A]A"="T>A", "G[T/A]C"="T>A", "G[T/A]G"="T>A", "G[T/A]T"="T>A", "T[T/A]A"="T>A", "T[T/A]C"="T>A", "T[T/A]G"="T>A",
	"T[T/A]T"="T>A", "A[T/C]A"="T>C", "A[T/C]C"="T>C", "A[T/C]G"="T>C", "A[T/C]T"="T>C", "C[T/C]A"="T>C", "C[T/C]C"="T>C", "C[T/C]G"="T>C", "C[T/C]T"="T>C", "G[T/C]A"="T>C", "G[T/C]C"="T>C",
	"G[T/C]G"="T>C", "G[T/C]T"="T>C", "T[T/C]A"="T>C", "T[T/C]C"="T>C", "T[T/C]G"="T>C", "T[T/C]T"="T>C", "A[T/G]A"="T>G", "A[T/G]C"="T>G", "A[T/G]G"="T>G", "A[T/G]T"="T>G", "C[T/G]A"="T>G",
	"C[T/G]C"="T>G", "C[T/G]G"="T>G", "C[T/G]T"="T>G", "G[T/G]A"="T>G", "G[T/G]C"="T>G", "G[T/G]G"="T>G", "G[T/G]T"="T>G", "T[T/G]A"="T>G", "T[T/G]C"="T>G", "T[T/G]G"="T>G", "T[T/G]T"="T>G"))

	#PLOT context re: Alexandrov
	ggp <- ggplot(df_ctxt, aes(x=context, y=norm, fill=color_ctxt))
	ggp + geom_histogram(stat="identity") + scale_fill_manual(values=c("#1EBBEB", "#000000", "#E0242A", "#A3A3A3", "#A3CC6E", "#EDC6C2")) + theme_bw() + theme(axis.text.x=element_text(angle = 90,vjust = 0.5))
}
