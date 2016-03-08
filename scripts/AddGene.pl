#!/usr/bin/perl
use strict;
use warnings;

open(ANN_FH, "$ARGV[0]");
my %GENE;
my %STRAND;
my $cols=0;
while(<ANN_FH>){
        chomp;
        my @local = split ("\t", $_);
	if (exists $GENE{"$local[0]\t$local[1]"}){
		$GENE{"$local[0]\t$local[1]"} = $GENE{"$local[0]\t$local[1]"}."; $local[3]";
		$STRAND{"$local[0]\t$local[1]"} = $STRAND{"$local[0]\t$local[1]"}."; $local[4]";
	}
	else{
		$GENE{"$local[0]\t$local[1]"} = $local[3];
		$STRAND{"$local[0]\t$local[1]"} = $local[4];
	}
}
close ANN_FH;

open (ORI,"$ARGV[1]");
while (<ORI>){
        chomp;
	my @temp = split("\t", $_);
	my $key  ="$temp[0]\t$temp[1]";
	if(exists $GENE{$key}){
		print "$_\t$GENE{$key}\t$STRAND{$key}\n";
	}
	else{
		print "$_\t--\t--\n";
	}
}
close ORI;
