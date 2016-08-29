#!/usr/bin/env Rscript
# This script is to make coveragePlot - this script is modfied from Rajesh's coverage.R

## example:
# transcript_coverage.R -f "/projects/Clinomics/ProcessedResults/OM161/20160415/CL0033_T3R_T/qc/CL0033_T3R_T.RnaSeqMetrics.txt /projects/Clinomics/ProcessedResults/OM161/20160415/CL0033_T_T/qc/CL0033_T_T.RnaSeqMetrics.txt /projects/Clinomics/ProcessedResults/OM161/20160415/CL0034_T_T/qc/CL0034_T_T.RnaSeqMetrics.txt " -o "test.png" -s "TEST"

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))
option_list <- list(
        make_option(c("-f", "--files"),  help="RnaSeqMetricsFile, required"),
        make_option(c("-o", "--image"),  help="output plot png file name, required"),
        make_option(c("-s", "--subject"),help="Title of the plot(subject name),required"),
        make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
        help="to output some information about a job.  [default: %default]")
    )
opt <- parse_args(OptionParser(option_list=option_list))

files  <- strsplit(opt$files, split="[\\s,;]", perl=TRUE)[[1]]
subject<-opt$subject

total <-length(files)+1

if (total > 8){
	cols <- brewer.pal(total, "Paired")
} else{
	cols <- brewer.pal(total, "Dark2")
}
png(opt$image, h=1200, w=1500, pointsize=25, type=c("cairo"))
#png(opt$image, h=1000, w=1000, pointsize=20)
par(mar = c(5, 5, 4, 13) + 0.1)
plot(c(-1,101), c(0,1.2), 
	type="n", 
	xaxt="n", 
	yaxt="n", 
	las=1,
	xlab='Normalized Distance Along Transcript', 
	cex.lab=1.5, 
	ylab='Normalized Coverage',
	main=subject
)
axis(1, cex.axis=1, at=c(0,20,40,60,80,100))
axis(2, cex.axis=1, at=c(0.0,0.2,0.4,0.6,0.8,1))
abline(  h=1, col="darkgray")
box()
labs<-c()
for (i in 1:length(files)){
	dat <- read.table(files[i], skip = 10, header =TRUE, sep ='\t')
	points(dat$normalized_position,dat$All_Reads.normalized_coverage,col=c(cols[i]))
	lines(dat$normalized_position,dat$All_Reads.normalized_coverage,col=c(cols[i]))
	#basename(files[i])
	labs<-c(labs,basename(files[i]))
}
par(xpd=TRUE)
labs <- paste("", gsub("Sample_|\\.RnaSeqMetrics.txt", "", labs, perl=TRUE), sep="")
#legend("topright", legend=labs, col=cols, lty=1, lwd=4, cex=0.5)
legend(106,1.2, legend=labs, col=cols, lty=1, lwd=4, cex=1)
