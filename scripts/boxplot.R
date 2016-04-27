library(stringr)
library(RColorBrewer)

args<-commandArgs(TRUE)
DIR = str_trim(args[1])
FILE=str_trim(args[2])
SAM=str_trim(args[3])

files <- list.files(path = DIR, pattern=".hotspot.depth$")

labs <- paste("", gsub("Sample_|\\.star|\\.bwa|\\.hotspot.depth", "", files, perl=TRUE), sep="")

cols <- brewer.pal(length(files)+1, "Dark2")

## y max
y_max =0 
for (i in 1:length(files)){
	cov.data   <-read.table(paste(DIR,files[i] ,sep = ""), sep="\t", quote="", head=F);
	cov.data <-cov.data[,6]
	if (max(cov.data) > y_max){ 
		y_max = max(cov.data)}
}
##

png(FILE,width = 1000, height = 1000, res=100, type=c("cairo"));
par(mar = c(18, 5, 4, 2) + 0.1)
plot(10,10,  xlim=c(1,length(files)+1), ylim=c(0,y_max+1), xaxt='n', yaxt="n", xlab='', ylab='Coverage', main=paste(SAM, "Hotspot Coverage", sep="\n"))

for (i in 1:length(files)){
        cov.data   <-read.table(paste(DIR,files[i] ,sep = ""), sep="\t", quote="", head=F);
	cov.data <-cov.data[,6]
	
	axis(1, labels=labs[i], at=i+0.5, las = 2, col=c(cols[i]), col.ticks =cols[i])
	boxplot(cov.data, add=T, at=i+0.5, las = 2, border=cols[i])
}
dev.off()
