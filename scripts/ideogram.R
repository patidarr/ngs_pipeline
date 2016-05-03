#!/usr/bin/env Rscript
# 
# ./ideogram.R -f Example.input -o test.pdf -t "Example"
# input contains 3 columns (CHR\tPOS\tP)
#
suppressPackageStartupMessages(library("pracma"))
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
	make_option(c("-f", "--file"), help="File containing list of numbers"),
	make_option(c("-o", "--out"), help="outputFile name(png/pdf)"),
	make_option(c("-t", "--title"), help="Main Title of the plot")
    )
opt <- parse_args(OptionParser(option_list=option_list))
file <-opt$f
if(grepl("png",opt$o)){
	png(opt$o, width = 1200, height = 500,type=c("cairo"))
}
if (grepl("pdf",opt$o)){
	pdf(opt$o, width=20,type=c("cairo"))
}
title=opt$t

ideogram <- function(x, chr="CHR", bp="POS", p="P", snp="SNP",
		col=c("gray10", "gray60"), logp=FALSE, ...) {
	# Check for sensible dataset
	## Make sure you have chr, bp and p columns.
	if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
	if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
	if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
	## make sure chr, bp, and p columns are numeric.
	if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
	if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
	if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
	# Create a new data.frame with columns called CHR, POS, and P.
	d=data.frame(CHR=x[[chr]], POS=x[[bp]], P=x[[p]])
	d <- subset(d, (is.numeric(CHR) & is.numeric(POS) & is.numeric(P)))
	d <- d[order(d$CHR, d$POS), ]
	#d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
	if (logp) {
		d$logp <- -log10(d$P)
	} else {
	d$logp <- d$P
	}
	d$pos=NA
	# Fixes the bug where one chromosome is missing by adding a sequential index column.
	d$index=NA
	ind = 0
	for (i in unique(d$CHR)){
		ind = ind + 1
		d[d$CHR==i,]$index = ind
	}
	nchr = length(unique(d$CHR))
	lastbase=0
	ticks=NULL
	xmal=NULL
	ablines=NULL	
	for (i in unique(d$index)) {
		if (i==1) {
			d[d$index==i, ]$pos=d[d$index==i, ]$POS
		} else {
			lastbase=lastbase+tail(subset(d,index==i-1)$POS, 1)
			d[d$index==i, ]$pos=d[d$index==i, ]$POS+lastbase
		}
		ticks = c(ticks, (min(d[d$CHR == i,]$pos) + max(d[d$CHR == i,]$pos))/2 + 1)
		xmax=max(d[d$CHR == i,]$pos)
		ablines=c(ablines, xmax)
	}
	xlabel = 'Human Genome (hg19)'
	ylabel="Log Relative Ratio"
	labs <- c(1:22,"X","Y")
	def_args <- list(xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=".",
		xlim=c(0,xmax), ylim=c(-2,+2),
		xlab=xlabel, ylab=ylabel)
	dotargs <- list(...)
	## And call the plot function passing NA, your ... arguments, and the default
	## arguments that were not defined in the ... arguments.
	do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
	# If manually specifying chromosome labels, ensure a character vector and number of labels matches number chrs.
	# Add an axis.
	axis(1, at=ticks, labels=labs, ...)
	# Create a vector of alternatiting colors
	col=rep(col, max(d$CHR))
	# Add points to the plot
	# if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
	icol=1
	for (i in unique(d$index)) {
		with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=".", ...))
		icol=icol+1
	}
	abline(v=ablines)
	abline(h=0)
	box()
}
data <- read.table(file, header=T, sep="\t")
ideogram(data,cex = 0.5, cex.axis = 0.8,col = c("blue4", "orange3"), main=title)
dev.off()
