#!/usr/bin/perl
use strict;
use warnings;

open(ANN_FH, "$ARGV[0]");
my %ANNOVAR;
my $cols=0;
while(<ANN_FH>){
        chomp;
        my @local = split ("\t", $_);
	if (exists $ANNOVAR{"$local[0]\t$local[1]"}){
		$ANNOVAR{"$local[0]\t$local[1]"} = $ANNOVAR{"$local[0]\t$local[1]"}."; $local[3]";
	}
	else{
		$ANNOVAR{"$local[0]\t$local[1]"} = $local[3];
	}
}
close ANN_FH;

open (ORI,"$ARGV[1]");
while (<ORI>){
        chomp;
	my @temp = split("\t", $_);
	my $key  ="$temp[0]\t$temp[1]";
	if(exists $ANNOVAR{$key}){
		print "$_\t$ANNOVAR{$key}\n";
	}
	else{
		print "$_\t--\n";
	}
}
close ORI;
