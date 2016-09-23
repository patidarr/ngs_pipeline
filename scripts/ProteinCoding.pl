#!/usr/bin/perl
use strict;
use warnings;
open(FH, $ARGV[0]);
while(<FH>){
	chomp;
	my @d  = split("\t",$_);
	if($_ =~ /^Chr/){
		print "$_\n";
	}
	elsif($d[5] =~ /^exonic/ and ($d[8] =~ /nonsynonymous/ or $d[8] =~ /stop/ or $d[8] =~ /frameshift/)){
		print "$_\n";
	}
	elsif($d[5] =~ /^splicing/){
		print "$_\n";
	}
}
