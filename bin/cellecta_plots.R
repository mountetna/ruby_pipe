#!/usr/bin/env Rscript

make_pval_plots = function(df_loc, sample_name, high_output_loc, low_output_loc){
    a = read.table(df_loc, sep='\t', header = T)
    write.table(a[order(a$gene),], df_loc, sep='\t', quote=F, row.names=F)
    a = read.table(df_loc, sep='\t', header = T)
    a$high_pval[a$high_pval == 0] = 1e-6
    a$low_pval[a$low_pval == 0] = 1e-6
    
    low = a[a$low_pval < 0.05,]
    low = low[order(low[,'low_pval'], decreasing = FALSE),]
    high = a[a$high_pval < 0.05,]
    high = high[order(high[,'high_pval'], decreasing = FALSE),]

    #plot high
    pdf(high_output_loc, height = 7, width = 7)
    par(mfrow=c(1,1),oma = c(0, 0, 3, 0))
    plot(a$high_pval, log='y', col= '#38B0DE50', cex = 0.7,
         main = 'Positive log fold change',
         xlab = "Gene(alphabetical)",
         ylab = "log(p-value)"         )
    text(row.names(high)[1:6], high$high_pval[1:6], labels = high$gene[1:6], cex = 0.6,srt = 35)
    mtext(sample_name, outer = TRUE, cex = 1.5)
    dev.off()

    #plot low
    pdf(low_output_loc, height = 7, width = 7)
    par(mfrow=c(1,1),oma = c(0, 0, 3, 0))
    plot(a$low_pval, log='y', col='#FFA50050', cex = 0.7,
         main = "Negative log fold change",
         xlab = "Gene(alphabetical)",
         ylab = "log(p-value)")
    text(row.names(low)[1:6], low$low_pval[1:6], labels = low$gene[1:6], cex = 0.6, srt = 35)

    mtext(sample_name, outer = TRUE, cex = 1.5)
    dev.off()

}

get_gene_list = function(gene_list_loc){
  con = file(gene_list_loc, open='r')
  genes = readLines(con)
  close(con)
  return(genes)
}

make_lfc_plots = function(df_loc, sample_name, output_loc, gene_list_loc){
    genes_of_interest = get_gene_list(gene_list_loc)
    print(paste('genes_of_interest', genes_of_interest))
    df = read.table(df_loc, sep='\t', header = T)
    luc_df = df[df$gene == 'Luc', ]

    pdf(output_loc, height = 7, width = 7)
        par( mfrow=c(1,1), oma = c(0,0,3,0) )
        palette('default')

        plot(log2(df$ctrl_count), log2(df$cond_count), cex = 0.8, col = '#99999970',
             xlab = 'log2(control counts)', ylab = 'log2(experimental counts)')
        points(log2(luc_df$ctrl_count), log2(luc_df$cond_count), col = 'black', cex = 1, pch = 20)
        
        #add points for genes of interest
        for(i in 1:length(genes_of_interest)){
            points(log2(df[df$gene == genes_of_interest[i], 'ctrl_count']),log2(df[df$gene == genes_of_interest[i], 'cond_count']),
            col = palette()[i+1],#i + 1 to skip black
            pch = 20, cex = 1 )
        }

        legend('topright', '', legend = c('luciferase', genes_of_interest), col = palette(), pch = 20, cex = 1)
        mtext(paste(sample_name, 'barcode counts'), outer = T, cex = 1.5)
    dev.off()

}

args = commandArgs(TRUE)
lib_dir = args[1] #i.e. the name of the R script file?
cmd = args[2] # i.e. name of function
print( as.list(args[-1:-2]))
do.call(cmd, as.list(args[-1:-2]))

