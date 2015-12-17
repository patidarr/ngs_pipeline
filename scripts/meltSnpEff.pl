#!/usr/bin/perl
use strict;
use warnings;
# This script parses the snpEff annotation for every transcript.
# Author Rajesh Patidar 
# 
# $0 input.snpeff.txt >output.snpeff.txt
# 
# caller --> order samples in vcf file --> snpEff --> formatInput --> This script --
#                                                                                    "
# caller --> order samples in vcf file --> snpEff --> formatInput --> This script ====>> Join Calls --> Annotation --> (extract columns of your interest)--> Filtering.
#                                                                                    "
# caller --> order samples in vcf file --> snpEff --> formatInput --> This script --


open (FH, $ARGV[0]);
my $Eff="Effect\tPutative_impact\tGene_Name\tGene_ID\tFeature_type\tFeature_ID\tTranscript_biotype\tHGVS.c\tHGVS.p";
while(<FH>){
	chomp;
	my @line  = split("\t", $_);
	my $partA = join("\t", @line[0..6]);
	my $partB = join("\t", @line[8..$#line]);
	if($_ =~ /^Chr/){
		print "$partA\tINFO\t$Eff\t$partB\n";
		next;
	}
	my $annotation  = $line[7];
	if($annotation =~ /;LOF=/){
		$annotation =~s/(.*);ANN\=(.*);LOF=(.*)//;
		my $info = $1;
		my $annotation = $2;
		my $lof = $3;
		if ($lof =~ /(.*);NMD=(.*)/){
			# LOF
			my @ann = split(",", $1);
			foreach my $sub (@ann){
				my ($Gene_Name, $Gene_ID, $Num_transcripts, $percent_affected) = processLOF($sub);
				print "$partA\t$info\tLoss-Of-Function\tHIGH\t$Gene_Name\t$Gene_ID\ttranscript\t$Num_transcripts\tNA\t\t$percent_affected\t$partB\n";
			}

			#NMD
			my @nmd = split(",", $2);
                        foreach my $element (@nmd){
                                my ($Gene_Name, $Gene_ID, $Num_transcripts, $percent_affected) = processLOF($element);
                                print "$partA\t$info\tNonsense-Mediated-Decay\tHIGH\t$Gene_Name\t$Gene_ID\ttranscript\t$Num_transcripts\tNA\t\t$percent_affected\t$partB\n";
                        }
		}
		else{
			# LOF
			my @ann = split(",", $lof);
			foreach my $sub (@ann){
				my ($Gene_Name, $Gene_ID, $Num_transcripts, $percent_affected) = processLOF($sub);
				print "$partA\t$info\tLoss-Of-Function\tHIGH\t$Gene_Name\t$Gene_ID\ttranscript\t$Num_transcripts\tNA\t\t$percent_affected\t$partB\n";
			}
		}	
		# Eff with LOF
		my @ann = split(",", $annotation);
		foreach my $sub (@ann){
			$sub = extractFields($sub);
			print "$partA\t$info\t$sub\t$partB\n";
		}	
	}
	else{
		# Eff
		$annotation =~s/(.*);ANN\=//;
		my $info = $1;	
		my @ann = split(",", $annotation);
		foreach my $sub (@ann){
			$sub = extractFields($sub);
			print "$partA\t$info\t$sub\t$partB\n";
		}
	}
}

sub extractFields{
	my ($input) = (@_);
	my ($Allele, $Effect, $Putative_impact, $Gene_Name, $Gene_ID, $Feature_type, $Feature_ID, $Transcript_biotype, $Rank, $HGVS_c, $HGVS_p, $cDNA_position, $CDS_position, $Protein_position, $Distance_to_feature, $informationMessage) = split /\|/, $input;
	return ("$Effect\t$Putative_impact\t$Gene_Name\t$Gene_ID\t$Feature_type\t$Feature_ID\t$Transcript_biotype\t$HGVS_c\t$HGVS_p");
}
sub processLOF{
	my ($input) = (@_);
	$input =~ s/\(//g;
	$input =~ s/\)//g;
	my ($Gene_Name, $Gene_ID, $Num_transcripts, $percent_affected) = split /\|/, $input;
	return ($Gene_Name, $Gene_ID, $Num_transcripts, $percent_affected);
}
