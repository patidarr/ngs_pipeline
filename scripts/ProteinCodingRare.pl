#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
local $SIG{__WARN__} = sub {my $message =shift; die $message;};

# Author: Rajesh Patidar rajbtpatidar@gmail.com
# This scripts take 4 inputs (in the order listed)
# 1	a list of variants (chr\tstart\tend\ref\talt\totherInfo) which should be filtered out from the 
# 2	a bed file of regions to be excluded from filtering out	
# 3	file to work on, must be annotated using khanlab annotation-pipeline/ngs-pipeline
# 4	MAF filter (0.05) This remove any variant if present in any manoj/minor  1000g/ES/ExAC population.
#
#
open(BLACK, $ARGV[0]);
open(WHITE, $ARGV[1]);
open(IN, $ARGV[2]);
my $freq = $ARGV[3];
my %REMOVE;
while(<BLACK>){
	chomp;
	my @t =split("\t", $_);
	$REMOVE{"$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]"} ="$t[5]";
}
close BLACK;
my %KEEP;
while(<WHITE>){
	chomp;
	my @t =split("\t", $_);
	$KEEP{"$t[0]\t$t[1]\t$t[2]"} ="$t[3]";
}
while(<IN>){
	chomp;
	if($_ =~ /^Chr/){
		print "$_\n"; 
		next;
	}
	my @a = split ("\t", $_);
	if (exists $REMOVE{"$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]"}){
		next;
	}
	foreach my $key (keys %KEEP){
		my @region= split("\t", $key);
		if($a[0] =~ /^$region[0]$/ and $a[1] >=$region[1] and $a[1] <= $region[2]){
			print "$_\n";
			next;
		}	
	}
	#if(exists $KEEP{"$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]"}){
	#	print "$_\n";
	#	next;
	#}
	if($a[5] =~ /^exonic/ and ($a[8] =~ /nonsynonymous/ or $a[8] =~ /stop/ or $a[8] =~ /frameshift/) or $a[5] =~ /^splicing/){
		for (@a[12..28]){
			s/^\.$/-2/;	
		}
		if( $a[12] <=$freq and $a[13] <=$freq and $a[14] <=$freq and $a[15] <=$freq and $a[16] <=$freq and $a[17] <=$freq and $a[18] <=$freq and $a[19] <=$freq and $a[20] <=$freq and $a[21] <=$freq and $a[22] <=$freq and $a[23] <=$freq and $a[24] <=$freq and $a[25] <=$freq and $a[26] <=$freq and $a[27] <=$freq and $a[28] <=$freq){
			print "$_\n";
		}
	}
}
close IN;
