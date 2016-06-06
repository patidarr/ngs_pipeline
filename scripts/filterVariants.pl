#!/usr/bin/perl
use strict;
use warnings;
my $freq = $ARGV[2];
open(IN, $ARGV[0]);
open(BLACK, $ARGV[1]);
my %LIST;
while(<BLACK>){
	chomp;
	my @t =split("\t", $_);
	$LIST{"$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]"} ="$t[5]";
}
while(<IN>){
	chomp;
	if($_ =~ /^Chr/){
		print "$_\n"; 
		next;
	}
	my @a = split ("\t", $_);
	if (exists $LIST{"$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]"}){
		next;
	}
	if( $a[12] <=$freq and $a[13] <=$freq and $a[14] <=$freq and $a[15] <=$freq and $a[16] <=$freq and $a[17] <=$freq and $a[18] <=$freq and $a[19] <=$freq and $a[20] <=$freq and $a[21] <=$freq and $a[22] <=$freq and $a[23] <=$freq and $a[24] <=$freq and $a[25] <=$freq and $a[26] <=$freq and $a[27] <=$freq and $a[28] <=$freq){
		print "$_\n";
	}
}
close IN;
