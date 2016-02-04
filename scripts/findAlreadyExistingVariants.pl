#!/usr/bin/perl
use strict;
use warnings;
## 
# Takes a list of query positions in file_1 and see if annotations are available in file_2 if they are write to file_3 else write to file_4
#
##
my %HASH;
open(ANN_FH, "$ARGV[0]");
while(<ANN_FH>){
        chomp;
        my @local = split ("\t", $_);
	$HASH{"$local[0]"}{$_} = '';
}
close ANN_FH;

open(BIG, "$ARGV[1]");
open(OUT1, ">$ARGV[2]");
while(<BIG>){
	chomp;
	my @local = split ("\t", $_);
	my $key = join "\t", @local[0..4];
	if(exists $HASH{"$local[0]"}{"$key"}){
		print OUT1 "$_\n";
		delete $HASH{"$local[0]"}{"$key"};
	}
}
close BIG;
close OUT1;
open(OUT2, ">$ARGV[3]");
print OUT2 "Chr\tStart\tEnd\tRef\tAlt\n";
foreach my $outter (keys %HASH) {
	foreach my $inner (keys %{$HASH{$outter}}) {
		print OUT2 "$inner\n";
	}
}
close OUT2;
