#!/opt/R/R-latest/bin/Rscript --no-save
suppressMessages(library(ABSOLUTE))

args = commandArgs(TRUE)

if (length(args) < 3) {
print("Usage: review_absolute.R <results_dir> <summary_name> <samples>")
quit()
}

results_dir = args[1]
summary_name = args[2]
samples = args[3]
samples = list.files(results_dir,samples)

CreateReviewObject(summary_name, file.path(results_dir,samples,paste(samples,".ABSOLUTE.RData",sep="")),
    file.path(results_dir, summary_name), "total")
