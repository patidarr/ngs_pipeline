library(stringr)
library(RColorBrewer)

args<-commandArgs(TRUE)
DIR = str_trim(args[1])
FILE=str_trim(args[2])
SAM=str_trim(args[3])

files <- list.files(path = DIR, pattern=".hotspot.depth$")

labs <- paste("", gsub("Sample_|\\.star|\\.bwa|\\.hotspot.depth", "", files, perl=TRUE), sep="")

cols <-c('#26294a','#01545a','#bd544f','#017351',
        '#03c383','#b8bd4f','#aad962','#fbbf45',
        '#bd8b4f','#ef6a32','#ed0346','#d76e60',
        '#a12a5e','#710162','#26294a','#01545a',
        '#bd544f','#017351','#03c383','#b8bd4f',
        '#aad962','#fbbf45','#bd8b4f','#ef6a32',
        '#ed0346','#d76e60','#a12a5e','#710162'
       )
png(FILE,width = 1000, height = 1000, res=100, type=c("cairo"));
par(mar = c(18, 5, 4, 2) + 0.1)
plot( c(1,length(files)+1), c(0,1000), type="n", xaxt="n", yaxt="n", las=1, xlab='', cex.lab=2, ylab='Coverage', main=paste(SAM, "Hotspot Coverage", sep="\n"))
axis(2, cex.axis=2, at=c(100,200,300,500,750,1000),las=2, cex.axis=1)
axis(2, cex.axis=2, at=c(50,400),las=2, cex.axis=1.2, col.axis="red")


for (i in 1:length(files)){
        cov.data   <-read.table(paste(DIR,files[i] ,sep = ""), sep="\t", quote="", head=F);
	cov.data <-cov.data[,6]
	axis(1, labels=labs[i], at=i+0.5, cex.axis=1.1,font=2,las = 2, col=c(cols[i]), col.ticks =cols[i])
	boxplot(cov.data, add=T, at=i+0.5, las = 2, border=cols[i],yaxt="n")
	med=print(round(median(cov.data)),digits=2)
	text(i+0.5,999,med,cex=1.5,font=2)
}
abline(h=c(50,400),col = "darkgray")
dev.off()
