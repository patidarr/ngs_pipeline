#!/usr/bin/env Rscript
# Takes as input a sample name to run sequenza for purity and poidy analysis
#

suppressPackageStartupMessages(library(sequenza))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
		make_option("--sample", help="Tumor sample name")
)
opt <- parse_args(OptionParser(option_list=option_list))
sample_name=opt$sample

DoSequenza <- function(sample_name) {
	message("run sequenza on sample = ", sample_name)

	chrom_list=c()
	for(i in c(1:22,'X','Y')){
			chrom_list=c(chrom_list, paste0('chr',i))
	}
	data_file <- paste(sample_name,".seqz_small.gz", sep="")
	data.file = data_file
	test <- sequenza.extract(data.file, chromosome.list=chrom_list)
	CP.example <- sequenza.fit(test)
	sequenza.results(sequenza.extract = test, cp.table = CP.example, sample.id = sample_name, out.dir=sample_name)
}
#sample_name= commandArgs(TRUE)[1]
log_file <- paste(sample_name,"_sequenza.log", sep="")
#sink(file = 'run_sequenza.log')
sink(file = log_file)
DoSequenza(sample_name)
sink()
