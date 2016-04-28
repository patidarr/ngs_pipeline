#!/bin/sh
#PBS -N ANNOVAR
set -eo pipefail
cd $1 ###
file="$2.anno"   ###
DATADIR=$4
CUSTOM=$3    ###
BUILD=hg19
###############################
# Add gene, cytoband,dbsnp, 1000g, ESP, CG69, NCI60 annotations
###############################
#       --dot2underline\
table_annovar.pl\
	$file\
	$DATADIR\
	-buildver ${BUILD}\
	-out $file\
	-remove\
	-protocol refGene,cytoBand,snp138,1000g2014oct_all,1000g2014oct_eur,1000g2014oct_afr,1000g2014oct_amr,1000g2014oct_eas,1000g2014oct_sas,esp6500_all,esp6500_ea,esp6500_aa,exac03nontcga,cg69,nci60\
	-operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f\
	-nastring "-1"
mv $file.hg19_multianno.txt $file.gene
rm -rf $file.refGene.invalid_input
sed -i '1s/\./_/g' $file.gene
###############################
# Add ExAC annotation
#
###############################
annotate_variation.pl\
	$file\
	$DATADIR\
	-buildver ${BUILD}\
	-otherinfo\
	-filter\
	-dbtype exac03

awk '{OFS="\t"};{print $3,$4,$5,$6,$7,$2}' $file.${BUILD}_exac03_dropped |sed -e 's/,/\t/g' >$file.exac.3
head -1 $DATADIR/${BUILD}_exac03.txt >>$file.exac.3
rm -rf $file.${BUILD}_exac03_dropped $file.${BUILD}_exac03_filtered
################################nno.sh
# Add clinseq annotation
#
################################
annotate_variation.pl\
	$file\
	$DATADIR\
	-buildver ${BUILD}\
	-otherinfo\
	-filter\
	-dbtype generic\
	-genericdbfile ${BUILD}_clinseq_951.txt

awk '{OFS="\t"};{print $3,$4,$5,$6,$7,$2}' $file.${BUILD}_generic_dropped |sed -e 's/,/\t/g' >$file.clinseq
head -1 $DATADIR/${BUILD}_clinseq_951.txt >>$file.clinseq
rm -rf $file.${BUILD}_generic_dropped $file.${BUILD}_generic_filtered
################################
# Add CADD annotation
#
################################
annotate_variation.pl\
	$file\
	$DATADIR\
	-buildver ${BUILD}\
	-otherinfo\
	-filter\
	-dbtype cadd


annotate_variation.pl\
	$file\
	$DATADIR\
	-buildver ${BUILD}\
	-otherinfo\
	-filter\
	-dbtype caddindel

cut -f 2-7 $file.${BUILD}_cadd_dropped $file.${BUILD}_caddindel_dropped |sed -e 's/,/\t/g' |awk '{OFS="\t"};{print $3,$4,$5,$6,$7,$1,$2}' >$file.cadd
head -1 $DATADIR/${BUILD}_caddindel.txt >>$file.cadd
rm -rf $file.${BUILD}_cadd_dropped $file.${BUILD}_cadd_filtered $file.${BUILD}_caddindel_dropped $file.${BUILD}_caddindel_filtered
################################
# Add Clinvar and COSMIC
#
################################
table_annovar.pl\
	$file\
	$DATADIR\
	-buildver ${BUILD}\
	--dot2underline\
	-out $file\
	-remove\
	-protocol cosmic76\
	-operation f\
	-nastring "NA"
mv $file.hg19_multianno.txt $file.cosmic
################################
# Add PCG 
#
################################
annotate_variation.pl\
	$file\
	$DATADIR\
	-buildver ${BUILD}\
	-otherinfo\
	-filter\
	-dbtype generic\
	-genericdbfile ${BUILD}_PCG_042616.txt

#	-genericdbfile ${BUILD}_PCG_112015.txt
awk -F "\t" '{OFS="\t"};{print $3,$4,$5,$6,$7,$2}' $file.${BUILD}_generic_dropped |sed -e 's/,/\t/g' >$file.pcg
head -1 $DATADIR/${BUILD}_PCG_042616.txt >>$file.pcg
rm -rf $file.${BUILD}_generic_dropped $file.${BUILD}_generic_filtered
################################
# Add HGMD
#
################################
OUT=`echo $file |sed -e 's/.anno//g'`
$CUSTOM $DATADIR/${BUILD}_clinvar_20160203.txt $file >$OUT.clinvar
################################
# Add HGMD
#
################################
OUT=`echo $file |sed -e 's/.anno//g'`
$CUSTOM $DATADIR/${BUILD}_hgmd.2016.1.txt $file >$OUT.hgmd
################################
# Add MATCH Trial
#
################################
$CUSTOM $DATADIR/${BUILD}_MATCHTrial_2015_11.txt $file >$OUT.match
################################
# Add MyCG
#${BUILD}_MCG.02.27.15.txt
#hg19_MCG.06.24.15.txt
################################
$CUSTOM $DATADIR/${BUILD}_MCG.06.24.15.txt $file >$OUT.mcg
################################
# Add DoCM
#
################################
$CUSTOM $DATADIR/${BUILD}_DoCM.txt $file >$OUT.docm
################################
################################
# Add CanDL
#
################################
$CUSTOM $DATADIR/${BUILD}_candl_10262015.txt $file >$OUT.candl
################################
# Add Targeted Cancer Care
#
################################
$CUSTOM $DATADIR/${BUILD}_targated_cancer_care_10262015.txt $file >$OUT.tcc
################################
# CiViC
#
################################
$CUSTOM $DATADIR/${BUILD}_civic_10262015.txt $file >$OUT.civic
################################
#
#
################################
rm -rf $file.invalid_input
rm -rf $file.refGene.invalid_input
rm -rf $file.log
rm -rf $file.anno.log 
