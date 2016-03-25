#!/usr/bin/env Rscript
# This script is calculate various stats on a given list (Only 1 list)
# This is written in mind for calculaing stats on a list coming from pipe "|" and stdout to file or pipe
## example: 
# ListStatistics.R -f SD -l file.input
# ListStatistics.R -f SD -l <(cat test)
#cat file.input| ListStatistics.R -f SD -l -
#cat file.input| ListStatistics.R -f SD -l /dev/stdin

#library("pracma")
suppressPackageStartupMessages(library("pracma"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

option_list <- list(
	make_option(c("-f", "--file"), help="File containing list of numbers"),
	make_option("--stat", help="Any of Mean, Median, SD, MAD, Min, Max, summary, sum")
    )
opt <- parse_args(OptionParser(option_list=option_list))

OpenRead <- function(arg) {

   if (arg %in% c("-", "/dev/stdin")) {
      file("stdin", open = "r")
   } else if (grepl("^/dev/fd/", arg)) {
      fifo(arg, open = "r")
   } else {
      file(arg, open = "r")
   }
}
file <-OpenRead(opt$f)
inData<-read.table(file, header=FALSE)



if(strcmp(toupper(opt$stat), 'SD') ) {
	out <-sd(inData$V1)
	cat(out, "\n")
}
if(strcmp(toupper(opt$stat), 'MAD')){
	out <-mad(inData$V1)
	cat(out, "\n")
}
if(strcmp(toupper(opt$stat), 'MEAN')){
	out <-mean(inData$V1)
	cat(out, "\n")
}
if(strcmp(toupper(opt$stat), 'AVERAGE')){
        out <-mean(inData$V1)
        cat(out, "\n")
}
if(strcmp(toupper(opt$stat), 'MEDIAN')){
        out <-median(inData$V1)
        cat(out, "\n")
}
if(strcmp(toupper(opt$stat), 'MIN')){
        out <-min(inData$V1)
        cat(out, "\n")
}
if(strcmp(toupper(opt$stat), 'MAX')){
        out <-max(inData$V1)
        cat(out, "\n")
}
if(strcmp(toupper(opt$stat), 'SUMMARY')){
        out <-summary(inData$V1)
	print(out)
#        cat(out, "\n")
}
if(strcmp(toupper(opt$stat), 'SUM')){
        out <-sum(inData$V1)
        cat(out, "\n")
}
