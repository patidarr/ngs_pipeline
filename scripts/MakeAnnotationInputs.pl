#!/usr/bin/perl
use warnings;
use strict;
my $file=$ARGV[0];

open(FH, "$file"); # Open User Input
open (ANFH, ">$file.anno"); # Write annovar input
open (FHpph, ">$file.pph");# write SIFT input
open (FHsift, ">$file.sift");# write PPH2 input 

while(<FH>){
	chomp;
	if($_ =~ /^#/ or $_ =~ /^Chr/ or $_ =~ /^CHR/ or $_ =~ /^chr\t/){
		print ANFH "Chr\tStart\tEnd\tRef\tAlt\n";
	}
	else{
		my @local = split ("\t", $_);
		if(defined $local[3] and $local[4]){
			print ANFH "$local[0]\t$local[1]\t$local[2]\t$local[3]\t$local[4]\n";
			if(length($local[3]) <2 and length($local[4]) <2 and  $local[3] =~ /[ATCG]{1}/ and $local[4] =~ /[ATCG]{1}/ and $local[3] !~ /-/ and $local[4] !~ /-/){
				print FHpph "$local[0]:$local[1]\t$local[3]/$local[4]\n";
				$local[0] =~ s/chr//;
				print FHsift "$local[0],$local[1],1,$local[3]/$local[4]\n"; #1,69094,1,G/C
			}
		}
	}
}
close FH;
close ANFH;
close FHpph;
close FHsift;
