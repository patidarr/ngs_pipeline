#!/usr/bin/perl
use strict;
use warnings;
open(FH, $ARGV[0]);
print "Chr\tStart\tEnd\tRef\tAlt\tGene\tHLA Allele\tPeptide Len\tSub-peptide Pos\tMutation Pos\tMT Epitope Seq\tWT Epitope Seq\tBest MT IC50 Method\tBest MT IC50\tCorresponding WT IC50\tCorresponding Fold Change\tMedian MT IC50\tMedian WT IC50\tMedian Fold Change\tNetMHC WT IC50\tNetMHC MT IC50\tNetMHCcons WT IC50\tNetMHCcons MT IC50\tNetMHCpan WT IC50\tNetMHCpan MT IC50\tPickPocket WT IC50\tPickPocket MT IC50\tSMM WT IC50\tSMM MT IC50\tSMMPMBEC WT IC50\tSMMPMBEC MT IC50\n";
while(<FH>){
	chomp;
	if($_ =~ /Chromosome/){
		next;
	}
	my @list =split("\t", $_);
	print "$list[0]\t$list[2]\t$list[2]\t$list[3]\t$list[4]\t";
	for (my $i=10; $i<=43; $i++){
		if (not defined $list[$i]){
			$list[$i] = "NA";
		}
		if ($list[$i] =~ /^\d+\.\d+/){
			$list[$i] =  sprintf("%.2f", $list[$i]);
		}
	}
	print join("\t",@list[10..20])."\t".join("\t",@list[29..43])."\n";
	#30-44
}
