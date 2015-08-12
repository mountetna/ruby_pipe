#Script to run babel in batch mode for use on computer cluster

args = commandArgs(TRUE)
lib_dir = args[1]
cmd = args[2]

if (length(args) < 1) {
  print("Usage:")
  print("  babel <lib_dir> doBabel <coverage_file> <output_file> <groups> <num_reps> <min_rna> <min_rpkm>")
  quit()
}

sample_names = function(groups,type=NA) {
	return(unlist(sapply(groups,function(group) {
		if (is.na(type)) {
			return( group$samples )
		}
		else if (type == "group") {
			return( rep( group$name, length(group$samples) ) )
		}
		else {
			return( paste(group$samples, type, sep=".") )
		}
	})))
}

run_babel = function(coverage,groups,min_rna,num_reps) {
	library(babel)

	names = c(
		sample_names(groups,"rna"),
		sample_names(groups,"rp")
	)

	babel.input = coverage[ , names ]

	table.RNA = babel.input[ , grep('.rna$',colnames(babel.input))]
	table.RP = babel.input[ , grep('.rp$',colnames(babel.input))]

	colnames(table.RNA) = sample_names(groups)
	colnames(table.RP) = sample_names(groups)

	babel.group = sample_names(groups,"group")

	set.seed(12345)

	babel.output = babel(table.RNA, table.RP, group=babel.group, nreps=num_reps,min.rna=min_rna)

	for (x in names(babel.output$within)) {
		babel.output$within[[x]]$Gene = rownames(coverage)
	}
	for (x in names(babel.output$between)) {
		babel.output$between[[x]]$Gene = rownames(coverage)
	}
	for (x in names(babel.output$combined)) {
		babel.output$combined[[x]]$Gene = rownames(coverage)
	}

	return(babel.output)
}

filter_by_rpkm = function(coverage,min_rpkm) {
	# this is a very rough way to get rpkm, don't sneeze too loud.

	gexp = coverage[ , c(-1,-2)]
	lib_size = colSums(gexp)
	rpkm = t(t( gexp / (coverage$size / 1e3) ) / (lib_size / 1e6))
	keep = apply(rpkm,1,function(gene) {
		return(any(gene > min_rpkm))
		})

	return(coverage[keep,])
}

make_groups = function(group_string) {
	groups = strsplit(group_string,":")[[1]]
	groups = lapply(groups, function(group) {
		group_list = strsplit(group,"=")[[1]]
		group_name = group_list[1]
		group_samples = strsplit(group_list[2],",")[[1]]
		return(list(name=group_name,samples=group_samples))
	})
	return(groups)
}

get_babel_coverage = function(coverage_file) {
	coverage = read.delim(coverage_file, header=TRUE, as.is=TRUE, sep="\t",row.names="gene_id")
	excluded =  c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique")
	coverage = coverage[ !rownames(coverage) %in% excluded, ]
	return(coverage)
}

doBabel = function(coverage_file, output_file, groups, num_reps,min_rna,min_rpkm) {
	coverage = get_babel_coverage(coverage_file)
	groups = make_groups(groups)
	num_reps = as.numeric(num_reps)
	min_rna = as.numeric(min_rna)
	min_rpkm = as.numeric(min_rpkm)
	filtered_coverage = filter_by_rpkm(coverage,min_rpkm)
	output.babel = run_babel(filtered_coverage, groups=groups, num_reps=num_reps,min_rna=min_rna)
	within.babel <- output.babel$within
	combined.babel <- output.babel$combined
	between.babel <- output.babel$between
	save("within.babel",file=paste(output_file,".within.babel", sep=""))
	save("combined.babel",file=paste(output_file,".combined.babel", sep=""))
	save("between.babel",file=paste(output_file,".between.babel", sep=""))
}

do.call(cmd,as.list(args[-1:-2]))
