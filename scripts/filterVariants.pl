#!/usr/bin/perl
use strict;
use warnings;
my $freq = $ARGV[1];
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	if($_ =~ /^Chr/){
		print "$_\n"; 
		next;
	}
	my @a = split ("\t", $_);
	if($a[12] <=$freq and $a[13] <=$freq and $a[14] <=$freq and $a[15] <=$freq and $a[16] <=$freq and $a[17] <=$freq and $a[18] <=$freq and $a[19] <=$freq and $a[20] <=$freq){
		print "$_\n";
	}
}
close IN;
