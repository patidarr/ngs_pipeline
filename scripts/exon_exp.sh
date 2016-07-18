#!/usr/bin/env bash
#
# Author: Rajesh Patidar rajbtpatidar@gmail.com
# Takes a bed file and calculate sum, RPKM and log2(RPKM) for every bed entry
#
#module load samtools/0.1.19
mil=1000000000
constant="1"
total_reads=$1
BedFile=$2
bamFile=$3
outFile=$4
if [ -f "$outFile" ];then
	rm $outFile
fi

#total_reads=`samtools flagstat $bamFile |head -1 | sed 's/\s/\t/g' | cut -f1`
#total_reads=407932018

while read -r line
do
	elements=( $line )
	chr=`echo ${elements[0]}|sed -e 's/chr//g'`
	mpileup_cord="$chr:${elements[1]}-${elements[2]}";
	sum=`samtools depth -Q 10 -r $mpileup_cord $bamFile | awk '{sum += $3} END {print sum}'`; 
	if [ -z "$sum" ];then
		sum="1"
		log=`echo "scale=3; l($sum)/l(2)" |bc -l`
		echo -e "$line\\t0\\t$sum\\t$log" >> $outFile # sum is 0
#		echo -e "$line\\t$sum" >> $outFile
	else
		RPKM=`echo "(($sum * $mil) / ($total_reads * (${elements[2]} - ${elements[1]}))) + $constant"|bc`
		log=`echo "scale=3; l($RPKM)/l(2)" |bc -l`
		echo -e "$line\\t$sum\\t$RPKM\\t$log" >> $outFile
#		echo -e "$line\\t$RPKM" >> $outFile
	fi
done < $BedFile
