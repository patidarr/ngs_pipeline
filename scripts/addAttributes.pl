#! /usr/local/bin/perl
use strict;
use warnings;
use 5.010;
local $SIG{__WARN__} = sub {
	my $message =shift;
	die $message;
};
my $hs      =`/bin/grep \"^$ARGV[0]\" $ARGV[1] |cut -f 21,22`;
chomp $hs;
my @cov =split("\t", $hs);
my $hs1  =sprintf("%.f", $cov[0]);
my $hs2 +=sprintf("%.f", $cov[1]);
my $in_file =$ARGV[2];
my $out_file=$ARGV[3];
open (IN, "<$in_file") or die "Can't open $in_file: $!\n";
open (OUT, ">$out_file") or die "Can't open $out_file: $!\n";
while(<IN>){
	chomp;
	my @line =split("\t", $_);
	if ($_ =~ /^$/){
		next;
	}
	if ($_ =~ /^#/){
		print OUT join("\t", @line[0..2]). "\tMEAN_BAIT_COVERAGE\tMEAN_TARGET_COVERAGE\t". join("\t", @line[3..$#line])."\n";
	}
	else{
		print OUT join("\t", @line[0..2]). "\t$hs1\t$hs2\t". join("\t", @line[3..$#line])."\n";
	}
}
close IN;
close OUT;
