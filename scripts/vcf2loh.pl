#!/usr/bin/perl
use strict;
use warnings;

open(FH , "$ARGV[0]");
print "Chr\tPosition\tLOH\n";
while(<FH>){
	chomp;
	if($_ =~ /DP4=(\d+),(\d+),(\d+),(\d+);/ and $_ !~ /INDEL/){
		if($1+$2+$3+$4 >=10){
			my @d = split("\t", $_);
			if($3+$4 eq "0"){
				print "$d[0]\t$d[1]\t0\n";	
			}
			else{
				my $vaf = sprintf("%2f",($3+$4)/($1+$2+$3+$4));
				print "$d[0]\t$d[1]\t$vaf\n";
			}
		}
	}
}
