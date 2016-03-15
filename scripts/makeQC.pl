#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dir="$ARGV[0]";
# Print Header
#

print "Sample\tLibrary\tDiagnosis\tTotal Reads\tRead Length\t%GC\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\tInsert Size (Median)\tMapping Quality (Mean)\tCoverage (Mean)\n";
#print "Sample\tTotal Reads\tRead Length\t%GC\tMapped Reads\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n";

for my $Sample( <$dir/*\/*>){
	if (-e "$Sample/qc"){
		my ($name,$path) = fileparse($Sample);
		my ($patient) =`basename $path`;
		chomp $patient;
		print "$patient\t";
		my $fastqc = "$path/$name/qc/fastqc/$name"."_R1_fastqc/fastqc_data.txt";	
		if (open (FH , $fastqc)){
			print "$name\t";
			my $diagnosis =`sed -n '/"Diagnosis"/, /}/ p' Diagnosis.json |grep "$name" |cut -d ":" -f 2 |sed -e 's/,//g'`;
			chomp $diagnosis;
			print "$diagnosis\t";
			while(<FH>){
				chomp;
				if($_ =~ /^Total Sequences\s+(\d+)/){
					print "$1\t";
				}
				if($_ =~ /Sequence length\s+(.*)/){
					print "$1\t";
				}
				if($_ =~ /%GC\s+(.*)/){
					print "$1\t";
				}
			}
			close FH;
		}
		my $picard = "$path/$name/qc/bwa.markdup.txt";
		if(open(FH, $picard)){
			while(<FH>){
				chomp;
				if ($_ =~ /^$name/){
					my @line  = split("\t", $_);
					shift(@line);	
					print join("\t", @line)."\t";
				}	
			}
			close FH;
		}
		elsif(open(FH, "$path/$name/qc/star.markdup.txt")){
			while(<FH>){
                                chomp;
                                if ($_ =~ /^$name/){
                                        my @line  = split("\t", $_);
                                        shift(@line);
                                        print join("\t", @line)."\t";
                                }
                        }
                        close FH;
		}
		my $bamqc = "$path/$name/qc/BamQC/genome_results.txt";
		if(open(FH, $bamqc)){
			while(<FH>){
				chomp;
				if($_ =~ /median insert size = (.*)/){
					print "$1\t";
				}
				if($_ =~ /mean mapping quality = (.*)/){
					print "$1\t";
				}
				if($_ =~ /mean coverageData = (.*)X/){
					print "$1"
				}
			}
		}
		print "\n";
	}
}	
