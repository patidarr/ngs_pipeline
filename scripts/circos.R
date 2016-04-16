suppressPackageStartupMessages(library(OmicCircos))
library(stringr)
library(RColorBrewer)

args<-commandArgs(TRUE)
DIR = str_trim(args[1])
FILE=str_trim(args[2])
SAM=str_trim(args[3])

files <- list.files(path = DIR, pattern=".loh$")

labs <- paste("", gsub("Sample_|\\.bwa|\\.star|\\.loh", "", files, perl=TRUE), sep="")

cols <- brewer.pal(length(files)+1, "Dark2")

#color<-c("salmon2","palevioletred1","rosybrown","red2","dodgerblue1","dodgerblue3");

options(stringsAsFactors = FALSE);
set.seed(1234);
png(FILE ,width = 1000, height = 1000, res=100, points=12, type=c("cairo"));
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

circos(R=400, cir="hg19", type="chr", mapping=UCSC.hg19.chr,print.chr.lab=TRUE, W=10, lwd=5, cex=1.5);

r=350
for (i in 1:length(files)){
        LOH.data   <-read.table(paste(DIR,files[i] ,sep = ""), sep="\t", quote="", head=T);
        circos(cir="hg19", R=r, W=50, type="s", mapping=LOH.data, col.v=3, col=cols[i], B=TRUE, cex=0.0001, lwd=1);
	r=r-45;
}
legend("topright", legend=labs, col=cols, lty=1, lwd=4, cex=0.5)
dev.off()
