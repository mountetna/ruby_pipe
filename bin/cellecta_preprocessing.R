#!/opt/R/R-latest/bin/Rscript --no-save

read.tsv = function(loc){
  return(read.table(loc, header = T, sep='\t', stringsAsFactors = F))
}

read_tables = function(ctrl_loc, cond_loc){
  control = read.tsv(ctrl_loc)
  colnames(control)[3] = 'ctrl_count'
  condition = read.tsv(cond_loc)
  colnames(condition)[3] = 'cond_count'
  # making the assumption that all barcodes are represented, 
  # and that tables are the same dimension.
  # otherwise will have to fill in missing bcs with reference and
  # will have to do something like this: merged[is.na(all_sets)] = 0
  merged = merge(control, condition, by= c('bc', 'gene'), all = T)
  merged = merged[!merged$gene == 'unknown',]
  merged$ctrl_count[is.na(merged$ctrl_count)] = 0
  merged$cond_count[is.na(merged$cond_count)] = 0
  return(merged)
}

normalize_data = function(df){
  # add 10 to thwart errors from 0
  df$ctrl_count = df$ctrl_count + 10
  df$cond_count = df$cond_count + 10
  df$ctrl_count = (df$ctrl_count/sum(df$ctrl_count))*1e6
  df$cond_count = (df$cond_count/sum(df$cond_count))*1e6

  return(df)
}

calc_log_fold_change = function(df){
  df$lfc = log2(df$cond/df$ctrl)
  return(df)
}

write_combined_table = function(df, out_file){
  write.table(df, out_file, sep='\t', quote = F, row.names = F)
}

normalize_and_log_fold_change = function(out_file, ctrl_loc, cond_loc){
  df = read_tables(ctrl_loc, cond_loc)
  df = normalize_data(df)
  df = calc_log_fold_change(df)
  write_combined_table(df, out_file)
}

args = commandArgs(TRUE)
lib_dir = args[1] #i.e. the name of the R script file?
cmd = args[2] # i.e. name of function

print( as.list(args[-1:-2]))


do.call(cmd, as.list(args[-1:-2]))