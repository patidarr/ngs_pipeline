#!/usr/bin/perl
use strict;
use warnings;
open(FH,"$ARGV[0]");
my %HASH;
shift @ARGV;
foreach my $file(@ARGV){
	if($file =~ /HC_DNASeq/){
		open(FH1,$file);
		while(<FH1>){
			chomp;
			my @a=split("\t", $_);
			my $key="$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[8]";
			if($a[7] =~ /FS=(.*?);/){
				$HASH{$key} =$1;
			}
		}
		close FH1;
	}
}

while(<FH>){
	chomp;
	my @a =split("\t", $_);
	if($a[8] =~ /HC_DNASeq/){
		my $key="$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]";
		$a[10] = $HASH{$key};
	}
	else{
		$a[10] = 'NA';
	}
	print join("\t",@a)."\n";
}
close FH;
