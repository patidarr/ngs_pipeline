#!/usr/bin/perl
use strict;
use warnings;

# Author Rajesh Patidar rajbtpatidar@gmail.com


my $file = $ARGV[0];
if(!$file){
	print STDERR "Give a vcf file\n";
	die;
}

unless (open(FH, "$file")){
	print STDERR "Can not open the give file $file\n";
	die;
}
while(<FH>){
	chomp;
	$_ =~ s/chr//g;
	my @local = split("\t", $_);
	if ($_ =~ /#/){ #Print column Header
		next;
	}
	my $depth =0;
        if ($local[7] =~ /DP=(\d+);/){
                $depth = $1;
        }
	if(length($local[3]) eq 1 and length($local[4]) eq 1 and $local[3] =~ /[ATCG]/ and $local[4] =~ /[ATCG]/ and $local[5] >=30 and $depth >=10){
		if ($local[9] =~ /0\/1/){
			print "$local[0]\t$local[1]\t$local[3]$local[4]\n";
		}
		elsif($local[9] =~ /1\/1/){
			print "$local[0]\t$local[1]\t$local[4]$local[4]\n";
		}
		elsif($local[9] =~ /0\/0/){
			print "$local[0]\t$local[1]\t$local[3]$local[3]\n";
		}
	}
}
close FH;
