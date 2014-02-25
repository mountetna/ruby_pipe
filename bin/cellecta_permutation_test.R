#!/opt/R/R-latest/bin/Rscript --no-save

# Find second highest log fold change for each gene

# take 10000 random samples of n barcodes and find second highest of that (where n is number of barcodes in your gene)
# however many times above number is higher than gene's calculated second best log fold change, divided by 10000 = p value

# do same with second lowest


order_hi_to_low = function(df){
  df = df[order(df[,5], decreasing = TRUE),]
  return(df)
}

get_second_best = function(df){
  return (df[2,5])
}
get_second_lowest = function(df){
  return(df[nrow(df)-1,5])
}

compare_rand_second_highs = function(df, second_high){
  sec_best_rand = get_second_best(df)
  if(sec_best_rand >= second_high){ return(1) }
  else{ return(0) }
}

compare_rand_second_lows = function(df, second_low){
  sec_low_rand = get_second_lowest(df)
  if(sec_low_rand <= second_low){ return(1) }
  else{ return(0) }
}


iterate = function(genes, curr_table){

  newlist = sapply(as.vector(genes), function(g) {
    set.seed(42234)

    print(paste("gene:", g, "start at", Sys.time()))

    curr_set = order_hi_to_low(curr_table[curr_table$gene == g, ]) #subset table for only one gene
    sec_best = get_second_best(curr_set)
    sec_worst = get_second_lowest(curr_set)

    count_higher = 0
    count_lower = 0
    num_permutations = 1:10000
    # num_permutations = 1:1000

    count_vec = as.data.frame(t(sapply(num_permutations, function(x){
      s = sample(nrow(curr_table), nrow(curr_set), replace=F)
      sub = curr_table[s,]
      sub = order_hi_to_low(sub)
      return( c( higher = compare_rand_second_highs(sub,sec_best), lower = compare_rand_second_lows(sub, sec_worst)) )
    })))

    pval_high = sum(count_vec$higher)/length(num_permutations)
    pval_low = sum(count_vec$lower)/length(num_permutations)
    print(paste("gene:", g, "complete at", Sys.time()))
    return( c( gene = g, second_highest = as.numeric(sec_best), high_pval = as.numeric(pval_high), second_lowest = sec_worst, low_pval = as.numeric(pval_low)) )
    }
   )
  df = as.data.frame(t(newlist))
  return(df)
}

read_gene_set = function(gene_set_loc){
  con = file(gene_set_loc, open='r')
  genes = readLines(con)
  close(con)
  return(genes)
}

permute = function(input_file, gene_set_loc, output_loc){
  df = read.table(input_file, sep='\t', header = TRUE)
  print(head(df))
  genes = read_gene_set(gene_set_loc)
  print(head(genes))
  print(paste('number of genes', length(genes)))
  results = iterate(genes, df)
  write.table(results, output_loc, sep='\t', quote=F, row.names=F)

}

args = commandArgs(TRUE)
lib_dir = args[1] #i.e. the name of the R script file?
cmd = args[2] # i.e. name of function
print( as.list(args[-1:-2]))
do.call(cmd, as.list(args[-1:-2]))




