#!/usr/bin/env Rscript

chr_str = function(v) {
	v = as.character(v)
	v[ v == "23" ] = "X"
	v[ v == "24" ] = "Y"
	paste("chr",v,sep="")
}

chr_int = function(v) {
	v = gsub("chr","",v)
	v[ v == "X" ] = "23"
	v[ v == "Y" ] = "24"
	as.integer(as.vector(v))
}

chr_sort=function(d) {
  return(d[order( chr_int(d[,1]) , d[,2]),])
}
