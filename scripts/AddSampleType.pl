#!/usr/bin/perl
use strict;
use warnings;
my $file = $ARGV[0];
my @sampleType = split(" ", $ARGV[1]);
print "\n\n\n";
my %HASH;

for (my $pair=0; $pair<=$#sampleType; $pair++){
	$HASH{$sampleType[$pair]} = $sampleType[$pair+1];
	$pair++;
}
print "$_\t$HASH{$_}\n" for (keys %HASH);

unless (open(IN, $file)){
	print STDERR "Can not open input file $file\n";
	die;
}
while(<IN>){
	chomp;
	my @line = split("\t", $_);
	print "$_\t".$HASH{$line[5]}."\n";
}
