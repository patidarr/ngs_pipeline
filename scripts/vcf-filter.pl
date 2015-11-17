#!/usr/bin/perl
use strict;
use warnings;
open(FH, $ARGV[0]);
while(<FH>){
	chomp;
	if($_ =~ /#/){
		print "$_\n";
	}
	else{
		my @line = split("\t", $_);
		if(grep {$_ eq '.' } @line[ 9 .. $#line ]){
#			print "$_\n"
		}
		else{
			print "$_\n";
		}	
	}
}
