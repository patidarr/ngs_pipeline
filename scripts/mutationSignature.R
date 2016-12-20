#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("deconstructSigs"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
                make_option("--input", help="input file name"),
                make_option("--sample", help="Header of the sample column"),
		make_option("--output", help="output pdf file name")
		
)
opt <- parse_args(OptionParser(option_list=option_list))
input=opt$input
output=opt$output
sample=opt$sample



mut_data <- read.table(input,sep="\t",header=T)

sigs.input <- mut.to.sigs.input(mut.ref = mut_data,
		sample.id = sample,
		chr = "Chr",
		pos = "Start",
		ref = "Ref",
		alt = "Alt")
sample_1 = whichSignatures(tumor.ref = sigs.input,
		signatures.ref = signatures.nature2013,
		sample.id = sample,
		contexts.needed = TRUE,
		tri.counts.method = 'default')

pdf(output)
plotSignatures(sample_1)
makePie(sample_1)
dev.off()
