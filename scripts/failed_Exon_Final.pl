#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
local $SIG{__WARN__} = sub { my $message =shift; die $message; };
#full path of coverage file
my $COV_FILE = $ARGV[0];
my $tmp = `basename $COV_FILE ".depth_per_base"`;
chomp $tmp;
my $OUTFILE = $ARGV[2];
my $OUTFILE2 = $ARGV[3];
my $THRES=$ARGV[1];
## feed coverage numbers into hash
my %cov_hash;
my %gene_hash;
my $prev_key;
my $base_cnt=0;
## open coverage file to store in hash
open(IN, $COV_FILE) or die "Couldn't open file $COV_FILE for processing: $!";
## skip first header line
<IN>;
while (<IN>) {
	chomp;
	my @cols = split(/\t/, $_);
	if ($cols[3] =~ /^RS\d+/){ # Skip Pharmacogenomics RSIDs available in Panel1 
		next;
	}
	my ($gene, $exnum) = split('___', $cols[3]);
	my $exon_key = join('_', $cols[0], $cols[1]);
	if (exists $cov_hash{$exon_key}) {
		$gene_hash{$gene}[0]++; # at gene length
		## check if depth below threshold, if yes increment count at [1] by 1 
		if($cols[5] < $THRES) { 
			$cov_hash{$exon_key}[1]++; # at exon level
			$gene_hash{$gene}[1]++; # at gene level
		}
	}
	else {
		## initialize hash with exon length and count=0 at [0] & [1] resp.
		my $exon_len = $cols[2] - $cols[1] + 1;
		$cov_hash{$exon_key}[0] = $exon_len;
		$cov_hash{$exon_key}[1] = 0;
		$cov_hash{$exon_key}[2] = $cols[3]; ## region name
		$cov_hash{$exon_key}[3] = $cols[0]; ## start
		$cov_hash{$exon_key}[4] = $cols[1]; ## end
		$cov_hash{$exon_key}[5] = $cols[2]; ## Chr
		## assign gene key and initialize the length
		$gene_hash{$gene}[0]++;
		$gene_hash{$gene}[1] = 0;
		if($cols[5] < $THRES) { 
			$cov_hash{$exon_key}[1]++; # at exon level
			$gene_hash{$gene}[1]++; # at gene level
		}
	}
}
close IN;

open(OUT, ">$OUTFILE") or die "Couldn't open file $OUTFILE for processing: $!";
print OUT "#Chr\tStart\tEnd\tName\tLength\t$tmp"."_".$THRES."x\n";
open(OUT2, ">$OUTFILE2") or die "Couldn't open file $OUTFILE2 for processing: $!";
print OUT2 "#Name\tLength\t$tmp"."_".$THRES."x\n";
## check each region
foreach my $exon_key (sort keys %cov_hash) {
	my $ratio = $cov_hash{$exon_key}[1]/$cov_hash{$exon_key}[0];
	if ( $ratio > 0.05 ) {
		my $line = join("\t", $cov_hash{$exon_key}[3], $cov_hash{$exon_key}[4], $cov_hash{$exon_key}[5], $cov_hash{$exon_key}[2], $cov_hash{$exon_key}[0], $ratio);
		print OUT "$line\n";
	}
} # end of foreach
foreach my $gene (sort keys %gene_hash) {
	my $ratio = $gene_hash{$gene}[1]/$gene_hash{$gene}[0];
	#print "$gene\t$gene_hash{$gene}[1]\t\t$gene_hash{$gene}[0]\t$ratio\n";
	if ( $ratio > 0.05 ) {
		my $line = join("\t", $gene, $gene_hash{$gene}[0], $ratio);
		print OUT2 "$line\n";
	}
} # end of foreach
close OUT;
close OUT2;
