#!/usr/bin/perl -ws
use List::Util qw(sum);
###########################
# Input is a file containing the window and snp file sorted using bedtools.
# module load bedtools
# awk '{OFS="\t"};{print $1,$2,$2,$3,$4}' $Sample >$Sample.bed
# bedtools makewindows -g /data/khanlab/ref/GATK/hg19/hg19.genome -w 100000 | awk '{OFS="\t"}{print $1,$2,$3,"win","end"}' >win.bed
# cat win.bed >>$Sample.bed
# sortBed -i $Sample.bed >$Sample.sorted.bed
# $0 $Sample.sorted.bed >$Sample.win.bed
#
#
#
###########################

#1       47851   47851   0.1085321       0.991861
#1       50251   50251   0.1211283       0
#1       51938   51938   0.326722        0
#1       52651   52651   0.08728445      0.002627225
#1       55338   55338   0.3408957       0.005440401
#1       64251   64251   0.1228615       0
$/="end\n";
open(FH, $ARGV[0]);
while(<FH>){
	chomp $_;
	my @line=split("\n", $_);
	pop @line;
	my @LRR;
	for (my $i=0; $i<=$#line; $i++){
		my @snp =split("\t", $line[$i]);
		push @LRR, $snp[3];
	}
	#my $mean =sprintf ("%.2f", mean(@LRR));
	my $mean =sprintf ("%.2f", Median(\@LRR));
	for (my $i=0; $i<=$#line; $i++){
		my @snp =split("\t", $line[$i]);
		print "$snp[0]\t$snp[1]\t$snp[2]\t$mean\n";
	}

}
sub mean{
	return @_ ? sum(@_) / @_ : 0
}
# usage Median(\@data);
##################################################
sub Median{
	my ($refdata) = @_;
	my $median;
	@$refdata = sort{$a<=>$b}@$refdata;
	my $count = @$refdata;
	if ($count <1){
		return (0);
	}
	if ($count %2){
		$median = $$refdata[int($count/2)];

	}
	else{
		$median = ($$refdata[$count/2]+ $$refdata[$count/2 -1])/2;
	}
#        $median =sprintf "%.6f", $median;
	return $median;
}
