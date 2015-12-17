#!/usr/bin/perl 
use strict;
use warnings;
#/data/khanlab/projects/patidar/Snakemake/scripts/MergeCalls.pl Sample_NCI0003tumor_E_C28PUACXX.MuTect.snpEff.txt Sample_NCI0003tumor_E_C28PUACXX.strelka.snvs.snpEff.txt Sample_NCI0003tumor_E_C28PUACXX.strelka.indels.snpEff.txt |./mergeMuTect_Strelka.pl - Sample_NCI0003tumor_E_C28PUACXX.MuTect.snpEff.txt |./mergeMuTect_Strelka.pl - Sample_NCI0003tumor_E_C28PUACXX.strelka.snvs.snpEff.txt >~/rrrr
#my $a = `head -1 @ARGV |sort |uniq |grep Chr |wc -l`;
#chomp $a;
#if ($a eq 1 ){
#	open (FH, "cut -f 1-5 @ARGV |sort |uniq |");
#}
#else{
#	print STDERR "order of samples in the input files does not match\n";
#	die;
#}
my %HoH;
for my $file (@ARGV){
	open(FH, $file);
	my $caller = getCaller($file);
	while(<FH>){
		chomp;
		my @line  = split("\t", $_);
		my $key = join "\t",@line[0..4];
		my $end = @line - 1 ;
		my $vcf =join "\t", @line[5..$end];
		$HoH{$key}{$caller} = $vcf;
	}
}

foreach my $snv (sort keys %HoH) {
	my $number = keys %{$HoH{$snv}};
	if ($snv =~ /Chr\tStart\tEnd/){
		print "$snv\t#Callers\tCallers\n";
		next;
	}
	print "$snv\t$number\t".join(",", (sort keys %{$HoH{$snv}}))."\n";
}
sub getCaller{
	my $caller = shift(@_);
	
	if ($caller =~ /platypus/i){
		return "Platypus";
	}
	elsif($caller =~ /mutect/i){
		return "MuTect";
	}
	elsif($caller =~ /Strelka/i){
		return("Strelka");
	}
	elsif($caller =~ /mpg/i){
		return ("Bam2MPG");
	}
	elsif($caller =~ /haplotype/i or $caller =~ /hapcaller/i){
		return ("HaplotypeCaller");
	}
	elsif($caller =~ /unified/i){
		return("UnifiedGenotyper");
	}
	elsif($caller =~ /freebayes/i){
		return ("FreeBayes");
	}
	else{
		return $caller;
	}
}

