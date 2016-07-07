#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
local $SIG{__WARN__} = sub {
	my $message =shift;
	die $message;
};
my $file = $ARGV[0];

########
my @sampleType = split(" ", $ARGV[1]);
my %HASH;
for (my $pair=0; $pair<=$#sampleType; $pair++){
	$HASH{$sampleType[$pair]} = $sampleType[$pair+1];
	$pair++;
}

########
my @capture = split(" ", $ARGV[2]);
my %HASH1;
for (my $pair=0; $pair<=$#capture; $pair++){
        $HASH1{$capture[$pair]} = $capture[$pair+1];
        $pair++;
}
########
unless (open(IN, $file)){
	print STDERR "Can not open input file $file\n";
	die;
}
while(<IN>){
	chomp;
	my @line = split("\t", $_);
	print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t".$HASH{$line[5]}."\t".$HASH1{$line[5]}."\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\n";
}
