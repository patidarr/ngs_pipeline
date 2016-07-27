#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
no if $] >= 5.017011, warnings => 'experimental::smartmatch';
use 5.010;
local $SIG{__WARN__} = sub {my $message =shift; die $message;};
open(FH, $ARGV[0]); # File Name
my $Patient =$ARGV[1];
my $Library =$ARGV[2];
my $Diagnosis=$ARGV[3];
#PF_BASES     PF_ALIGNED_BASES  RIBOSOMAL_BASES  CODING_BASES  UTR_BASES  INTRONIC_BASES  INTERGENIC_BASES  IGNORED_READS  CORRECT_STRAND_READS  INCORRECT_STRAND_READS  PCT_RIBOSOMAL_BASES  PCT_CODING_BASES  PCT_UTR_BASES  PCT_INTRONIC_BASES  PCT_INTERGENIC_BASES  PCT_MRNA_BASES  PCT_USABLE_BASES  PCT_CORRECT_STRAND_READS  MEDIAN_CV_COVERAGE  MEDIAN_5PRIME_BIAS  MEDIAN_3PRIME_BIAS  MEDIAN_5PRIME_TO_3PRIME_BIAS  SAMPLE  LIBRARY  READ_GROUP
#10472336775  10471765851       35636850         19039801      2419240    1370225         10413301463       0              0                     0                       0.003403             0.001818          0.000231       0.000131            0.994417              0.002049        0.002049          0   
print "#Patient\tLibrary\tDiagnosis\tALIGNED_BASES\tRIBOSOMAL_BASES\tCODING_BASES\tUTR_BASES\tINTRONIC_BASES\tINTERGENIC_BASES\tPCT_RIBOSOMAL_BASES\tPCT_CODING_BASES\tPCT_UTR_BASES\tPCT_INTRONIC_BASES\tPCT_INTERGENIC_BASES\tPCT_MRNA_BASES\tPCT_USABLE_BASES\n";
while(<FH>){
	chomp;
	next if $_=~ /^$/ or $_ =~ /^#/ or $. >9;
	my @a = split("\t", $_);
	if ($#a >3 and $a[0] =~ /\d+/){
		print "$Patient\t$Library\t$Diagnosis\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]";
		foreach (@a[10..16]){
			$_ = $_*100 if $_ =~ /\d+/;
			print "\t$_";
		}
		print "\n";
	}
}
close FH;
