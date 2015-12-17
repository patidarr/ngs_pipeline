#!/usr/bin/env Rscript
# This script is to make coveragePlot - this script is modfied from Rajesh's coverage.R

## example: 
# coveragePlot.R -f "CL0023_N_P.final.bam.depth CL0024_N_P.final.bam.depth" -o "test"
# opt=NULL
# opt$covFiles="SUBJECT/32/32N/ucsc.hg19.bwamem/qc/32N.final.bam.depth SUBJECT/32/32T/ucsc.hg19.bwamem/qc/32T.final.bam.depth"
# opt$outName = 'SUBJECT/32/ucsc.hg19.bwamem/32_clin.ex.v1.coveragePlot.png'
# coveragePlot.R -f "SUBJECT/32/32N/ucsc.hg19.bwamem/qc/32N.final.bam.depth SUBJECT/32/32T/ucsc.hg19.bwamem/qc/32T.final.bam.depth" -o  SUBJECT/32/ucsc.hg19.bwamem/32_clin.ex.v1.coveragePlot.png

suppressPackageStartupMessages(library("optparse"))
option_list <- list( 
	make_option(c("-f", "--covFiles"), help="histFile, required"),
	make_option(c("-o", "--outName"), default='coveragePlot.png', help="output name. [default: %default]"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
	help="to output some information about a job.  [default: %default]")		
    )

# get command line options, if help option encountered print help and exit,
opt <- parse_args(OptionParser(option_list=option_list))

if( ! is.element('covFiles', names(opt)) ) stop("Option for covFiles is required. ")

files <- strsplit(opt$covFiles, split="[\\s,;]", perl=TRUE)[[1]]

cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.delim(paste(files[i], sep=""), header=FALSE, sep="\t", as.is=TRUE)
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][, 5])
}
labs <- basename( gsub("\\.final\\.bam\\.depth", "", files, perl=TRUE) )

library(RColorBrewer)
cols <- brewer.pal(length(cov)+1, "Dark2")

# Save the graph to a file
png(opt$outName, h=1000, w=1000, pointsize=20)
# Create plot area, but do not plot anything. Add gridlines and axis labels.
par(las=2)
plot(
	c(1,1000),
	c(0,1.001),
	log="x",
	type='n', 
	xaxt="n",
	yaxt="n",
	xlab="Coverage", 
	ylab="Fraction of capture target bases \u2265 depth", 
	ylim=c(0,1.001),
	las=1,
	xlim=c(1,1000),
	main="Target Region Coverage"
	)
axis(1, at=c(5,10,20,50, 100, 200, 500, 1000))
axis(2, at=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
abline(  v=c(5,10,20,50, 100, 200, 500, 1000), col = "darkgray")
abline(  h=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), col="darkgray")

# Actually plot the data for each of the alignments (stored in the lists).
for (i in 1:length(cov)){
	points(
		cov[[i]][1:1000, 2], cov_cumul[[i]][1:1000], type='l', lwd=3, col=cols[i]
	)
}
box()
# Add a legend using the nice sample labeles rather than the full filenames.
legend("topright", legend=labs, col=cols, lty=1, lwd=4)
dev.off()
