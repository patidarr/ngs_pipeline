#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
local $SIG{__WARN__} = sub {
        my $message =shift;
        die $message;
};
#########################
# Author: Rajesh Patidar rajbtpatidar@gmail.com
#
# This script takes a formatted vcf file (tdf) (Somatic) and calls on RNASeq of the same tumor and try to find
#  if the somatic indels are present in the RNASeq or not (Expressed).
#
#########################
my $SOMATIC = $ARGV[0];
my $RNA     = $ARGV[1];
my $RNAname = $ARGV[2];
open (IN1, $SOMATIC);
open(IN2, $RNA);


my %Expressed;
while(<IN2>){
	chomp;
	my @line = split("\t",$_);
	$Expressed{join("\t",@line[0..4])} = join("\t",@line[9..13]);
}

while(<IN1>){
	chomp;
	if ($_ =~ /^Chr/){
		print "$_\t$RNAname.GT\tRNASeq.TotCov\tRNASeq.RefCov\tRNASeq.VarCov\tRNASeq.VAF\n";
	}
	else{
		my @line = split("\t",$_);
		if (exists $Expressed{join("\t",@line[0..4])}){
			print "$_\t".$Expressed{join("\t",@line[0..4])}."\n";
		}
		else{
			my ($status,$info) = search($line[0], $line[1], $line[2]);
			if ($status eq 'match'){
				print "$_\t$info\n";
			}
			else{
				print "$_\t0\t0\t0\t0\t0\n";
					
			}
		}
	}
}
sub search{
	my ($chr, $start, $end) = (@_);
	foreach my $indel (keys %Expressed){
		my @ind = split("\t", $indel);
		if ($ind[0] eq $chr and ($start ~~ [$ind[1] - 10 .. $ind[1] + 10]  or $end ~~ [$ind[2] - 10 .. $ind[2] + 10])){
			return("match", $Expressed{$indel});
		}
	}
}
