#!/usr/bin/perl
use strict;
use warnings;
open(ANN_FH, "$ARGV[0]");
my %ANNOVAR;
while(<ANN_FH>){
	chomp;
	my @local = split ("\t", $_);
	my $end = @local - 3 ;
	my $key = join "\t", @local[0..$end];
	if ($_ =~ /^Chr/){
                print "#$key\tLevel_Somatic\tLevel_Germline\n";
                next;
        }
	if(not exists $ANNOVAR{"$key"}){
		$ANNOVAR{"$key"} = $local[-1];
	}
}
close ANN_FH;

open (ORI,"$ARGV[1]");
while (<ORI>){
	chomp;
	my @local = split ("\t", $_);
	my $end = @local - 3 ;
	my $key = join "\t", @local[0..$end];
	if ($_ =~ '^Chr'){
		next;
	}
	if (exists $ANNOVAR{$key}){
		print "$key\t$local[-1]\t$ANNOVAR{$key}\n";
		delete $ANNOVAR{$key};
	}
	else{
		 print "$key\t$local[-1]\tNA\n";
	}
}
foreach (keys %ANNOVAR){
	print "$_\tNA\t$ANNOVAR{$_}\n";
}
