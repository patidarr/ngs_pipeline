#!/bin/sh
set -eo pipefail

# Author: Rajesh Patidar rajbtpatidar
# 	This script add Gene, Transcript, Exon strand info to a 3 column bed file. 
# 	The output is compatible with the pipeline this reside in.
# 	If your input do not contain "chr" in first column please add it before running me.
#	Also remove any header and unwanted chromosomes. contigs..



dir=`dirname $0`

BED=$1
OUT=$2
REFSEQ=$dir/ref/Exons_06302016_refseq_formatted_noUTR.merged.bed
ENS=$dir/ref/Exons_06302016_ensembl_formatted_noUTR.merged.bed

module load bedtools
module load igvtools


# Copy the input as design file, should be replaced by the vendor provided design file whereever available
cp $BED ${OUT}.refseq.ensembl.noUTR.notfound.design.hg19.bed
igvtools index ${OUT}.refseq.ensembl.noUTR.notfound.design.hg19.bed



slopBed -i $BED -g $dir/ref/ucsc.hg19.genome -b 100 >${OUT}
intersectBed -a $REFSEQ -b ${OUT} -wa >${OUT}.foundInRefSeq
intersectBed -a ${OUT} -b $REFSEQ -v >${OUT}.PaddedNotfoundInRefSeq
intersectBed -a $BED -b ${OUT}.PaddedNotfoundInRefSeq -wa -u > ${OUT}.NotfoundInRefSeq
intersectBed -a $ENS -b ${OUT}.NotfoundInRefSeq >${OUT}.foundInENSEMBL #remove -wa option
intersectBed -a ${OUT}.NotfoundInRefSeq -b $ENS -v >${OUT}.NotfoundInENSEMBL


awk '{OFS="\t"}{print $1,$2,$3,"NOTFOUND",0,"NULL"}' ${OUT}.NotfoundInENSEMBL >${OUT}.Notfound
# remove the padding from remaining regions, mergeBed and make Gene and Exon numbers ___ seperated.
cat ${OUT}.foundInRefSeq ${OUT}.foundInENSEMBL ${OUT}.Notfound |sortBed -i - |uniq|mergeBed -i - -c 4,5,6 -o distinct,collapse,collapse |awk '{OFS="\t"}{print $1,$2,$3,$4"___"$5,$6}' >${OUT}.refseq.ensembl.noUTR.notfound.target.hg19.bed
igvtools index ${OUT}.refseq.ensembl.noUTR.notfound.target.hg19.bed

slopBed -i ${OUT}.refseq.ensembl.noUTR.notfound.target.hg19.bed -g $dir/ref/ucsc.hg19.genome -b 20 >${OUT}.refseq.ensembl.noUTR.notfound.targetbp.hg19.bed
rm -rf ${OUT} ${OUT}.foundInRefSeq ${OUT}.NotfoundInRefSeq ${OUT}.foundInENSEMBL ${OUT}.NotfoundInENSEMBL ${OUT}.Notfound ${OUT}.PaddedNotfoundInRefSeq
igvtools index ${OUT}.refseq.ensembl.noUTR.notfound.targetbp.hg19.bed
rm -rf igv.log
