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
cosmic = whichSignatures(tumor.ref = sigs.input,
		signatures.ref = signatures.cosmic,
		sample.id = sample,
		contexts.needed = TRUE,
		tri.counts.method = 'default')

nature = whichSignatures(tumor.ref = sigs.input,
                signatures.ref = signatures.nature2013,
                sample.id = sample,
                contexts.needed = TRUE,
                tri.counts.method = 'default')

pdf(output,width=10)
plotSignatures(cosmic, sub='Mutational Signature Based on COSMIC')
makePie(cosmic, sub='Mutational Signature Based on COSMIC')
plotSignatures(nature, sub='Mutational Signature based on Nature 2013--23945592')
makePie(nature, sub='Mutational Signature based on Nature 2013--23945592')

dev.off()
