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
                make_option("--resultOut", help="outputFile for Normalized FPKM")
)
opt <- parse_args(OptionParser(option_list=option_list))
threads=opt$nt
libName=opt$lib
target_file=opt$targetFile	
referenceGTF_file=opt$referenceGTF
featureType=opt$featureType
annotationRDS = opt$annotationRDS
resultOut=paste(opt$resultOut, "/", sep="")


fpkmToTpm <- function(fpkm){
    exp(log(fpkm) - log(sum(fpkm,na.rm =TRUE)) + log(1e6))
}


if(toupper(featureType) == "GENE") GTFAttrType="gene_id"
if(toupper(featureType) == "TRANSCRIPT") GTFAttrType="transcript_id"
if(toupper(featureType) == "EXON") GTFAttrType="exon_id"

# count numbers of reads mapped to reference genome at Gene Transcript & Exon
fc <- featureCounts(files=target_file,annot.ext=referenceGTF_file,isGTFAnnotationFile=TRUE,GTF.featureType="exon" ,GTF.attrType=GTFAttrType,useMetaFeatures=TRUE,allowMultiOverlap=FALSE,nthreads=threads,isPairedEnd=TRUE,requireBothEndsMapped=FALSE,countChimericFragments=TRUE)
colnames(fc$counts) <- libName
colnames(fc$stat) <-libName
##Get count & annotation Obejct
countObj=fc$counts
Annotation=readRDS(annotationRDS) %>% data.frame()

########################################### Choose the Annotation
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
GeneDF.CPM   <- as.data.frame(cpm(GeneDF.Norm,  normalized.lib.sizes = TRUE,log = FALSE))  ; colnames(GeneDF.CPM) <- c(libName)
GeneDF.rpkm  <- as.data.frame(rpkm(GeneDF.Norm, normalized.lib.sizes = TRUE, log = FALSE)) ; colnames(GeneDF.rpkm) <- c(libName)
GeneDF.tpm   <- apply(GeneDF.rpkm, 2 , fpkmToTpm)       ; colnames(GeneDF.tpm) <- c(libName)

########################################### Prepare final files
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
filePrefix = paste(resultOut,libName, ".", featureType, sep="")
saveRDS(fc, file=paste(filePrefix, "fc.RDS",sep="."))
write.table(GeneDF_Norm_CPM, file=paste(filePrefix, "CPM.txt",sep="."), sep="\t",row.names = FALSE, quote = FALSE)
write.table(GeneDF_Norm_rpkm, file=paste(filePrefix, "FPKM.txt",sep="."), sep="\t",row.names = FALSE, quote = FALSE)
write.table(GeneDF_Norm_tpm, file=paste(filePrefix, "TPM.txt",sep="."), sep="\t",row.names = FALSE, quote = FALSE)
