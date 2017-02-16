library(stringr)
args<-commandArgs(TRUE)
DIR = str_trim(args[1])
FILE=str_trim(args[2])
SAM=str_trim(args[3])


files <- list.files(path = DIR, pattern="coverage.txt$")

cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(paste(DIR,"/",files[i], sep=""))
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
}
labs <- paste("", gsub("Sample_|\\.bwa|\\.star|\\.coverage\\.txt", "", files, perl=TRUE), sep="")

library(RColorBrewer)

cols <-c('#26294a','#01545a','#bd544f','#017351',
        '#03c383','#b8bd4f','#aad962','#fbbf45',
        '#bd8b4f','#ef6a32','#ed0346','#d76e60',
        '#a12a5e','#710162','#26294a','#01545a',
        '#bd544f','#017351','#03c383','#b8bd4f',
        '#aad962','#fbbf45','#bd8b4f','#ef6a32',
        '#ed0346','#d76e60','#a12a5e','#710162'
       )
# Save the graph to a file
png(FILE, h=1000, w=1000, pointsize=10, type=c("cairo"))

# Create plot area, but do not plot anything. Add gridlines and axis labels.
par(las=2,cex.lab=2.2, cex.main=2.4, mar=c(6,7,6,2), mgp=c(4,1,0))
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
	cex=3,
	xlim=c(1,1000),
	main=paste(SAM, "\nTarget Region Coverage", sep="\t")
	)
axis(1, cex.axis=2, at=c(5,10,20,50, 100, 200, 500, 1000))
axis(2, cex.axis=2, at=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
abline(  v=c(5,10,20,50, 100, 200, 500, 1000), col = "darkgray")
abline(  h=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), col="darkgray")

# Actually plot the data for each of the alignments (stored in the lists).
for (i in 1:length(cov)){
	points(
		cov[[i]][2:1001, 2], cov_cumul[[i]][1:1000], type='l', lwd=3, col=cols[i]
	)
}
box()
# Add a legend using the nice sample labeles rather than the full filenames.
legend("bottomleft", legend=labs, col=cols, lty=1, lwd=4, cex=2)

dev.off()
