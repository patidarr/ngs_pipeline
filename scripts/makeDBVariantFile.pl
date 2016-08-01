#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use List::Util qw(first);
use 5.010;
local $SIG{__WARN__} = sub {
	my $message =shift;
	die $message;
};
#
# Author Rajesh Patidar rajbtpatidar@gmail.com
# This script takes the Patient ID and make a table of all the variants in the following format to be loaded to the database
# Chr\tStart\tEnd\tRef\tAlt\tCaller\tQUAL\tFS\tSample\tTotalReads\tAltReads    ## For every library not patient.
#  Caller, QUAL and FS are comma seperated
#
# Will take a list of file in following order.
#
# muTect expressed or mutect
# Strelka Expressed or Strelka (SNPS)
# Strelka Indels
# Haplotype Caller 
# Platypus
# MPG
# Haplotype Caller RNASeq
# If a variant is found in first 3 files it should be rejected from next 3 files. (Somatic)
# ./makeDBVariantFile.pl NCI0002/Sample_NCI0002tumor_E_C4PJ0ANXX/calls/Sample_NCI0002tumor_E_C4PJ0ANXX.MuTect.annotated.expressed.txt NCI0002/Sample_NCI0002tumor_E_C4PJ0ANXX/calls/Sample_NCI0002tumor_E_C4PJ0ANXX.strelka.snvs.annotated.expressed.txt NCI0002/Sample_NCI0002tumor_E_C4PJ0ANXX/calls/Sample_NCI0002tumor_E_C4PJ0ANXX.strelka.indels.annotated.txt NCI0002/NCI0002/calls/NCI0002.hapcaller.annotated.txt NCI0002/NCI0002/calls/NCI0002.platypus.annotated.txt NCI0002/NCI0002/calls/NCI0002.bam2mpg.annotated.txt NCI0002/Sample_NCI0002tumor_T_D1UTYACXX/calls/Sample_NCI0002tumor_T_D1UTYACXX.hapcaller.snpEff.txt
#
#my $last_index_of_annotation = 4;
my $last_index_of_annotation = 157;
my $index_qual = 158;
my $index_fs   = 158;
my (%CALLER, %QUAL, %FS, %READS);

foreach my $file(@ARGV){
	unless (open(FH, $file)){
		print STDERR "Can not open the file $file\n";
		exit;
	}
	my @a = split(/[.]/, basename($file));
	my $caller = $a[1];
	my %samples;
	while(<FH>){
		chomp;
		my @line = split("\t", $_);
		if($. ==1 and $_ =~ /^Chr/){
			$last_index_of_annotation = first { $line[$_] eq 'ACMG_LSDB' } 0..$#line;
			$index_qual = first { $line[$_] eq 'QUAL' } 0..$#line;
			$index_fs   = first { $line[$_] eq 'QUAL' } 0..$#line;	
			for (my $i=$last_index_of_annotation+5; $i<$#line; $i+=5){
				my $tmp = $line[$i];
				$tmp =~ s/\.GT//;
				$samples{$i} = $tmp;
			}
			next;
		}
		if($line[0] =~ /_/){
			next;
		}
		my $key  = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]";
		foreach (keys %samples){
			if ( not exists $READS{"$key\t$samples{$_}"}){
				$READS{"$key\t$samples{$_}"} = "$line[$_+1]\t$line[$_+3]";
				$QUAL{"$key\t$samples{$_}"}  = "$line[$index_qual]";
				$FS{"$key\t$samples{$_}"}  = "$line[$index_fs]";
				$CALLER{"$key\t$samples{$_}"}  = "$caller";
				
			}
			else{
				$QUAL{"$key\t$samples{$_}"} = $QUAL{"$key\t$samples{$_}"}.";$line[$index_qual]";
				$FS{"$key\t$samples{$_}"} = $FS{"$key\t$samples{$_}"}.";$line[$index_fs]";
				$CALLER{"$key\t$samples{$_}"} = $CALLER{"$key\t$samples{$_}"}.";$caller";
			}
		}
	}
	close FH;
	#foreach (keys %samples) { print "$_\t$samples{$_}\n";}	
}
# Print the data.
foreach (sort keys %READS) {
	my ($total, $alt) =split("\t", $READS{$_});
	if ($total =~ /^\d+$/ and $alt =~ /^\d+$/){
		print "$_\t$CALLER{$_}\t$QUAL{$_}\t$FS{$_}\t$READS{$_}\n";
	}
}
