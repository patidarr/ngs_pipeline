#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
local $SIG{__WARN__} = sub {
        my $message =shift;
        die $message;
};

unless (open(FH, $ARGV[0])){
	print STDERR "Can not find file $ARGV[0]";
	exit;
}
print "Chr\tStart\tEnd\tRef\tAlt\tSIFT Prediction\tSIFT Score\n";
while(<FH>){
	chomp;
	my @data = split("\t", $_);
	if ($. == 1){
	}
	else{
		$data[0] =~ s/\//,/g;
		my @head = split(",", $data[0]);	
		if($data[10] =~ /EXON/){
			print "chr$head[0]\t$head[1]\t$head[1]\t$head[3]\t$head[4]\t$data[13]\t$data[14]\n";
		}
	}
}
close FH;
