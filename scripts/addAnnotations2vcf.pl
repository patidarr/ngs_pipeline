#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
local $SIG{__WARN__} = sub {
        my $message =shift;
        die $message;
};
open(ANN_FH, "$ARGV[0]");
my %ANNOVAR;
while(<ANN_FH>){
	chomp;
	my @local = split ("\t", $_);
	my $key = join "\t", @local[0..4];
	my $end = @local - 1 ;
	my $value = join "\t", @local[5..$end];
	if(not exists $ANNOVAR{"$key"}){
		$ANNOVAR{"$key"} = $value;
	}
}
close ANN_FH;

open (ORI,"$ARGV[1]");

while (<ORI>){
	chomp;
	my @temp = split("\t", $_);
	my $val;
	my $vcf;	
	my $end = @temp - 1 ;
	$val = "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]";
	$vcf = join "\t", @temp[5..$end];
	if (exists $ANNOVAR{$val}){
		print "$val\t$ANNOVAR{$val}\t$vcf\n";
	}
}
