#!/usr/bin/perl
use strict;
use warnings;

open(ANNOTATION, "$ARGV[0]");
my %HASH;
my $cols=0;
while(<ANNOTATION>){
        chomp;
	$_ =~ s/^#//g;
        my @local = split ("\t", $_);
	my $key = join "\t", @local[0..4];
	my $end = @local - 1 ;
	$cols = $end - 5 + 1;
	my $value = join "\t", @local[5..$end];
	if( not exists $HASH{"$key"}){
		$HASH{"$key"} = $value;
	}
}
close ANNOTATION;

open (FILE,"$ARGV[1]");
my $string="\t-1"x$cols;
while (<FILE>){
        chomp;
	my @temp = split("\t", $_);
	my $val =join "\t", @temp[0..4];
	my $end = @temp - 1 ;
	my $vcf =join "\t", @temp[5..$end];
	if (exists $HASH{$val}){
		print "$_\t$HASH{$val}\n";
	}
	else{
		print "$_"."$string\n";
	}
		
}
close FILE
