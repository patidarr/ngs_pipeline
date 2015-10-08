#!/bin/sh
#PBS -N anno-pipe

CODE=$1
DIR=$2
FILE=$3
DATADIR="/data/Clinomics/Ref/annovar/"
TOOL="$CODE/addAnnotation.pl"
cd $DIR
$TOOL $FILE.anno.exac.3   $FILE.anno.gene >$FILE.anno.gene.exac
$TOOL $FILE.anno.clinseq  $FILE.anno.gene.exac >$FILE.anno.gene.exac.clinseq
$TOOL $FILE.anno.cadd     $FILE.anno.gene.exac.clinseq >$FILE.anno.gene.exac.clinseq.cadd
$TOOL $FILE.sift.out      $FILE.anno.gene.exac.clinseq.cadd >$FILE.anno.gene.exac.clinseq.cadd.sift
$TOOL $FILE.polyphen2.out $FILE.anno.gene.exac.clinseq.cadd.sift >$FILE.anno.gene.exac.clinseq.cadd.sift.pph
$TOOL $FILE.anno.clinvar  $FILE.anno.gene.exac.clinseq.cadd.sift.pph >$FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar
$TOOL $FILE.hgmd     $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar >$FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd
$TOOL $FILE.match    $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd >$FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match
$TOOL $FILE.docm     $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match >$FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm 
$TOOL $FILE.mcg      $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm >$FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg
$TOOL $FILE.anno.pcg           $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg >$FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg.pcg
$TOOL $FILE.uvm      $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg.pcg >$FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg.pcg.uvm
$TOOL $FILE.germline      $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg.pcg.uvm >$FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg.pcg.uvm.germline
$CODE/GeneAnnotation.pl $DATADIR/hg19_ACMG.txt $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg.pcg.uvm.germline >$FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg.pcg.uvm.germline.acmg

cp $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg.pcg.uvm.germline.acmg $FILE.annotations.txt

$CODE/joinAllAnnotations.pl $FILE.annotations.txt $FILE >$FILE.annotations.final.txt


#cat SIFT.e* SIFT.o* PPH2.e* PPH2.o* ANNOVAR.e* ANNOVAR.o* swb*.o* swb*.e* pph_roundup.e* pph_roundup.o*
rm -rf SIFT.e* SIFT.o* PPH2.e* PPH2.o* ANNOVAR.e* ANNOVAR.o* swb*.o* swb*.e* pph_roundup.e* pph_roundup.o*
rm -rf pph_roundup_*.e pph_roundup_*.o pph_swarm_*.e pph_swarm_*.o
rm -rf $FILE.anno $FILE.pph $FILE.sift $FILE.sift_predictions.tsv $FILE.sift.out $FILE.pph2 $FILE.anno.gene $FILE.anno.exac.3 $FILE.anno.clinseq  $FILE.pph2.out $FILE.polyphen2.out $FILE.anno.cadd $FILE.anno.clinvar $FILE.anno.pcg $FILE.hgmd $FILE.match $FILE.mcg $FILE.docm $FILE.germline $FILE.uvm

rm -rf $FILE.anno.gene.exac $FILE.anno.gene.exac.clinseq $FILE.anno.gene.exac.clinseq.cadd $FILE.anno.gene.exac.clinseq.cadd.sift $FILE.anno.gene.exac.clinseq.cadd.sift.pph $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg.pcg $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg.pcg.uvm $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg.pcg.uvm.germline $FILE.anno.gene.exac.clinseq.cadd.sift.pph.clinvar.hgmd.match.docm.mcg.pcg.uvm.germline.acmg 

echo -e "Annotation Pipeline finished on File $FILE\n Please collect result from $DIR \n\nRegards, \nOncogenomics Section\nNCI" |mutt  -s "Annotation Pipeline Status"  `whoami`@mail.nih.gov -c patidarr@mail.nih.gov
