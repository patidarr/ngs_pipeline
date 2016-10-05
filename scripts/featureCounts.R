#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("Rsubread"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("pracma"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
		make_option("--nt", help="Number of threads"),
		make_option("--lib", help="Library name to be printed in output files"),
                make_option("--targetFile", help="path of Bam"),
                make_option("--referenceGTF", help="GTF file"),
                make_option("--countOut", help="outputFile for read Count"),
                make_option("--fpkmOut", help="outputFile for Normalized FPKM")
)
opt <- parse_args(OptionParser(option_list=option_list))
threads=opt$nt
library=opt$lib
target_file=opt$targetFile
referenceGTF_file=opt$referenceGTF
count_file=opt$countOut
fpkm_file=opt$fpkmOut

# count numbers of reads mapped to reference genome at Gene Transcript & Exon
fc_Gene <- featureCounts(files=target_file,annot.ext=referenceGTF_file,isGTFAnnotationFile=TRUE,GTF.featureType= "exon",GTF.attrType="gene_id",useMetaFeatures=TRUE,allowMultiOverlap=FALSE,nthreads=threads,isPairedEnd=TRUE,requireBothEndsMapped=FALSE,countChimericFragments=TRUE,reportReads=FALSE)
fc_Transcript <- featureCounts(files=target_file,annot.ext=referenceGTF_file,isGTFAnnotationFile=TRUE,GTF.featureType= "exon",GTF.attrType="transcript_id",useMetaFeatures=TRUE,allowMultiOverlap=TRUE,nthreads=threads,isPairedEnd=TRUE,requireBothEndsMapped=FALSE,countChimericFragments=TRUE,reportReads=FALSE)
fc_Exon <- featureCounts(files=target_file,annot.ext=referenceGTF_file,isGTFAnnotationFile=TRUE,GTF.featureType= "exon",GTF.attrType="exon_id",useMetaFeatures=TRUE,allowMultiOverlap=TRUE,nthreads=threads,isPairedEnd=TRUE,requireBothEndsMapped=FALSE,countChimericFragments=TRUE,reportReads=FALSE)

#Save Binary count files of Gene Transcript & Exon
saveRDS(fc_Gene,file=paste(count_file,".Gene.fc.RDS",sep=""))
saveRDS(fc_Transcript,file=paste(count_file,".Transcript.fc.RDS",sep=""))
saveRDS(fc_Exon,file=paste(count_file,".Exon.fc.RDS",sep=""))

#Make DF
xGene <- DGEList(counts=fc_Gene$counts, genes=fc_Gene$annotation[,c("GeneID","Length")])
xTranscript <- DGEList(counts=fc_Transcript$counts, genes=fc_Transcript$annotation[,c("GeneID","Length")])
xExon <- DGEList(counts=fc_Exon$counts, genes=fc_Exon$annotation[,c("GeneID","Length")])

## Normalize and FPKM
x_rpkm_preGene <- calcNormFactors(xGene) ; x_rpkm_preTranscript <- calcNormFactors(xTranscript) ; x_rpkm_preExon <- calcNormFactors(xExon)
x_rpkm_Gene <- rpkm(x_rpkm_preGene,x_rpkm_preGene$genes$Length) ; x_rpkm_Transcript <- rpkm(x_rpkm_preTranscript,x_rpkm_preTranscript$genes$Length) ; x_rpkm_Exon <- rpkm(x_rpkm_preExon,x_rpkm_preExon$genes$Length) 

##Get Annotation
annotation_dfGene <- data.frame(t(apply(fc_Gene$annotation, 1, function(x) { chr <- as.character(unlist(strsplit(x[2],";",fixed=TRUE)))[1] ;start <- min(as.numeric(unlist(strsplit(x[3], ";", fixed=TRUE))));gene <-x[1] ;end <- max(as.numeric(unlist(strsplit(x[4], ";", fixed=TRUE))));return(c(chr,start,end,gene))})), stringsAsFactors=FALSE)
annotation_dfTranscript <- data.frame(t(apply(fc_Transcript$annotation, 1, function(x) { chr <- as.character(unlist(strsplit(x[2],";",fixed=TRUE)))[1] ;start <- min(as.numeric(unlist(strsplit(x[3], ";", fixed=TRUE))));gene <-x[1] ;end <- max(as.numeric(unlist(strsplit(x[4], ";", fixed=TRUE))));return(c(chr,start,end,gene))})), stringsAsFactors=FALSE)
annotation_dfExon <- data.frame(t(apply(fc_Exon$annotation, 1, function(x) { chr <- as.character(unlist(strsplit(x[2],";",fixed=TRUE)))[1] ;start <- min(as.numeric(unlist(strsplit(x[3], ";", fixed=TRUE))));gene <-x[1] ;end <- max(as.numeric(unlist(strsplit(x[4], ";", fixed=TRUE))));return(c(chr,start,end,gene))})), stringsAsFactors=FALSE)

##Reformat dfs
count_dfGene = cbind(annotation_dfGene,data.frame(fc_Gene$annotation[,c("Length")],fc_Gene$counts,stringsAsFactors=FALSE)) 
count_dfTranscript = cbind(annotation_dfTranscript,data.frame(fc_Transcript$annotation[,c("Length")],fc_Transcript$counts,stringsAsFactors=FALSE)) 
count_dfExon = cbind(annotation_dfExon,data.frame(fc_Exon$annotation[,c("Length")],fc_Exon$counts,stringsAsFactors=FALSE)) 
colnames(count_dfGene) <-  c("Chr", "Start", "End" , "GeneID", "Length", library)
colnames(count_dfTranscript) <-  c("Chr", "Start", "End" , "TranscriptID", "Length", library)
colnames(count_dfExon) <-  c("Chr", "Start", "End" , "ExonID", "Length", library)

fpkm_dfGene = cbind(annotation_dfGene,data.frame(x_rpkm_Gene,stringsAsFactors=FALSE))
fpkm_dfTranscript = cbind(annotation_dfTranscript,data.frame(x_rpkm_Transcript,stringsAsFactors=FALSE))
fpkm_dfExon = cbind(annotation_dfExon,data.frame(x_rpkm_Exon,stringsAsFactors=FALSE))
colnames(fpkm_dfGene) <-  c("Chr", "Start", "End" , "GeneID", library)
colnames(fpkm_dfTranscript) <-  c("Chr", "Start", "End" , "TranscriptID", library)
colnames(fpkm_dfExon) <-  c("Chr", "Start", "End" , "ExonID", library)

# write rpkm and count file
write.table(count_dfGene,file=paste(count_file,".Gene.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
write.table(count_dfTranscript,file=paste(count_file,".Transcript.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
write.table(count_dfExon,file=paste(count_file,".Exon.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
write.table(fpkm_dfGene,file=paste(fpkm_file,".Gene.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
write.table(fpkm_dfTranscript,file=paste(fpkm_file,".Transcript.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
write.table(fpkm_dfExon,file=paste(fpkm_file,".Exon.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
