#!/usr/bin/perl
use strict;
use warnings;
#Gene_refGene    ACMG_Disease    ACMG_Age-to-Report      ACMG_Gene-Reviews-PubMedID      ACMG_Inheritance        ACMG_Known-vs-Expected  ACMG_LSDB
open(ANNOTATION, "$ARGV[0]");
my %HASH;
my $cols=0;
while(<ANNOTATION>){
        chomp;
	$_ =~ s/^#//g;
        my @local = split ("\t", $_);
	$cols = @local;
	$HASH{"$local[0]"} = "$_";
}
close ANNOTATION;
open (FILE,"$ARGV[1]");
my $string="\t-"x$cols;
while (<FILE>){
        chomp;
	my @temp = split("\t", $_);
	my $val =$temp[$ARGV[2]];
	if (exists $HASH{$val}){
		print "$_\t$HASH{$val}\n";
	}
	else{
		print "$_"."$string\n";
	}
		
}
close FILE
