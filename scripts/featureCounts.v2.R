#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("Rsubread"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("pracma"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
                make_option("--nt", help="Number of threads"),
                make_option("--lib", help="Library name to be printed in output files"),
                make_option("--targetFile", help="path of Bam"),
                make_option("--referenceGTF", help="GTF file"),
		make_option("--featureType", help="gene or transcript or exon"),
		make_option("--annotationRDS", help="UCSC or ENSEMBL RDS "),
                make_option("--countOut", help="outputFile for read Count"),
                make_option("--resultOut", help="outputFile for Normalized FPKM")
)
opt <- parse_args(OptionParser(option_list=option_list))
threads=opt$nt
library=opt$lib
target_file=opt$targetFile	
referenceGTF_file=opt$referenceGTF
featureType=opt$featureType
annotationRDS = opt$annotationRDS
count_file=opt$countOut
exprssion_file=opt$resultOut


fpkmToTpm <- function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


if(toupper(featureType) == "GENE") GTFAttrType="gene_id"
if(toupper(featureType) == "TRANSCRIPT") GTFAttrType="transcript_id"
if(toupper(featureType) == "EXON") GTFAttrType="exon_id"

print(paste(GTFAttrType,annotationRDS))
# count numbers of reads mapped to reference genome at Gene Transcript & Exon
fc <- featureCounts(files=target_file,annot.ext=referenceGTF_file,isGTFAnnotationFile=TRUE,GTF.featureType="exon" ,GTF.attrType=GTFAttrType,useMetaFeatures=TRUE,allowMultiOverlap=FALSE,nthreads=threads,isPairedEnd=TRUE,requireBothEndsMapped=FALSE,countChimericFragments=TRUE,reportReads=FALSE)

#fc <- readRDS("/data/khanlab/apps/featureCountTPM.FPKM.CPM/counts.UCSC.transcript.fc.RDS")

#Save Binary count files of Gene Transcript & Exon
saveRDS(fc, file=paste(count_file,annotationType, featureType, "fc.RDS",sep="."))

##Get count & annotation Obejct
countObj=fc$counts
Annotation=readRDS(annotationRDS)

########################################### Choose the Annotation
print(paste("Chooseing the Annotation"))
if( toupper(featureType) == "EXON" )
{
  rownames(Annotation) <- Annotation$ExonID
  Annotation$ExonID    <- factor(Annotation$ExonID, levels=rownames(countObj))
  Annotation           <- Annotation %>% dplyr::arrange(ExonID)
  
  genesObj <- Annotation[,c("ExonID", "Length")]
} else if( toupper(featureType) == "TRANSCRIPT") {
  rownames(Annotation) <- Annotation$TranscriptID
  Annotation$TranscriptID <- factor(Annotation$TranscriptID, levels=rownames(countObj))
  Annotation           <- Annotation %>% dplyr::arrange(TranscriptID)
  
  genesObj <- Annotation[,c("TranscriptID", "Length")]
} else{
  rownames(Annotation) <- Annotation$GeneID
  Annotation$GeneID <- factor(Annotation$GeneID, levels=rownames(countObj))
  Annotation           <- Annotation %>% dplyr::arrange(GeneID)
  
  genesObj <- Annotation[,c("GeneID", "Length")]
}

########################################### Choose the Method 
##  Make EdgeR Object
colnames(genesObj) <- c("GeneID", "Length")
GeneDF_EdgeR       <- DGEList(counts=countObj, genes=genesObj)
## Estimate Normalising Factors
GeneDF.Norm  <- calcNormFactors(GeneDF_EdgeR) ; 
## Regularized Log Transformation using CPM, FPKM & TPM values
GeneDF.CPM   <- as.data.frame(cpm(GeneDF.Norm,  normalized.lib.sizes = TRUE,log = FALSE))
GeneDF.rpkm  <- as.data.frame(rpkm(GeneDF.Norm, normalized.lib.sizes = TRUE, log = FALSE))
GeneDF.tpm   <- apply(rpkm(GeneDF.Norm, normalized.lib.sizes = TRUE), 2 , fpkmToTpm)

########################################### Prepare final files
print(paste("Prepare final files"))
if( toupper(featureType) == "EXON" )
{
  GeneDF_Norm_CPM  <- cbind(data.frame(Annotation[,c("Chr","Start","End","GeneName","TranscriptID","ExonID")]), GeneDF.CPM )
  GeneDF_Norm_rpkm <- cbind(data.frame(Annotation[,c("Chr","Start","End","GeneName","TranscriptID","ExonID")]), GeneDF.rpkm )
  GeneDF_Norm_tpm  <- cbind(data.frame(Annotation[,c("Chr","Start","End","GeneName","TranscriptID","ExonID")]), GeneDF.tpm )
  
} else if( toupper(featureType) == "TRANSCRIPT" ){
  
  GeneDF_Norm_CPM  <- cbind(data.frame(Annotation[,c("Chr","Start","End","GeneName","TranscriptID")]), GeneDF.CPM )
  GeneDF_Norm_rpkm <- cbind(data.frame(Annotation[,c("Chr","Start","End","GeneName","TranscriptID")]), GeneDF.rpkm )
  GeneDF_Norm_tpm  <- cbind(data.frame(Annotation[,c("Chr","Start","End","GeneName","TranscriptID")]), GeneDF.tpm )
  
} else{
  GeneDF_Norm_CPM  <- cbind(data.frame(Annotation[,c("Chr","Start","End","GeneName")]), GeneDF.CPM)
  GeneDF_Norm_rpkm <- cbind(data.frame(Annotation[,c("Chr","Start","End","GeneName")]), GeneDF.rpkm)
  GeneDF_Norm_tpm  <- cbind(data.frame(Annotation[,c("Chr","Start","End","GeneName")]), GeneDF.tpm)
}

########################################### Choose approprite folder and write files
write.table(GeneDF_Norm_CPM, file=paste(exprssion_file, featureType, "CPM.txt",sep="."), sep="\t",row.names = FALSE, quote = FALSE)
write.table(GeneDF_Norm_rpkm, file=paste(exprssion_file, featureType, "FPKM.txt",sep="."), sep="\t",row.names = FALSE, quote = FALSE)
write.table(GeneDF_Norm_tpm, file=paste(exprssion_file, featureType, "TPM.txt",sep="."), sep="\t",row.names = FALSE, quote = FALSE)
