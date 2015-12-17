#!/bin/env Rscript
## This script is to order sample columns in a vcf file 
## Author: Jack Zhu@2015-12-08 

## example: 
# vcfOrderCol.R -i NCI0201.bam2mpg.raw.vcf -o NCI0201.bam2mpg.raw.ordered.vcf
# unset R_HOME to get rid of the following 
# WARNING: ignoring environment value of R_HOME

suppressPackageStartupMessages(library("optparse"))
option_list <- list( 
	make_option(c("-i", "--inVCF"), help="input vcf file, required"),
	make_option(c("-o", "--outVCF"), default='', help="output sample column ordered vcf. [default: %default]"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
	help="to output some information about a job.  [default: %default]")		
    )

opt <- parse_args(OptionParser(option_list=option_list))
if( ! is.element('inVCF', names(opt)) ) stop("Option for inVCF is required. ")

inVCF <- opt$inVCF
outVCF <- opt$outVCF
if( outVCF == '' ) {
    outVCF = sub('.vcf', '.ordered.vcf', inVCF)
}

df <- read.delim( pipe( paste('grep -v "^##" ', inVCF, sep="") ), sep="\t", as.is=TRUE, comment.char = "", check.names=FALSE)

headers <- names(df)
formatCol <- which(headers == "FORMAT")
if( ncol(df) > formatCol ) {
    headersOrdered <- c(headers[1:formatCol], sort(headers[ (formatCol+1): ncol(df) ]))
    df <- df[ , headersOrdered ]
} 

options(warn=-1)
system( paste('grep "^##" ', inVCF, " > ", outVCF, sep="") )
write.table(df, file=outVCF, append=TRUE, sep="\t", na="", row.names=FALSE, quote=FALSE)

