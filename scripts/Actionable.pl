#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(first);

if($ARGV[0] eq 'somatic'){
	Somatic($ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4]);
	# 1 == HotspotFile [reference]
	# 2 == CancerGene Census [reference]
	# 3 == somaticFile
	# 4 == Annotations
}
elsif($ARGV[0] eq 'germline' or $ARGV[0] eq 'rnaseq'){
	Germline($ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4], $ARGV[5]);
	# 1 == somaticFile
	# 2 == germlineFile
	# 3 == Annotations
	# 4 == CancerGene List [reference]
	# 5 == HotspotFile [reference]
}
else{
	die $!;
}
sub Germline{
	# First remove anything somatic
	my ($somatic, $germline, $annotation, $cancerGeneList, $hotspot) = (@_);

	unless (open(ANN_FH, "$somatic")){
		print STDERR "Can not open file $somatic\n";
		exit;
	}
	# Somatic file to hash
	my %DATA;
	while(<ANN_FH>){
		chomp;
		my @local = split ("\t", $_);
		my $key = join "\t", @local[0..4];
		$DATA{"$key"} = "somatic";
	}
	close ANN_FH;
	# Annotations to hash
	unless (open(REF, "$annotation")){
		print STDERR "Can not open file $annotation\n";
		exit;
	}
	my %ANNOTATION;
	while(<REF>){
		chomp;
		my @local = split ("\t", $_);
		my $key = join "\t", @local[0..4];
		my $end = @local - 1 ;
		my $value = join "\t", @local[5..$end];
		$ANNOTATION{$key} =$value;
	}
	close REF;
	# CancerGeneList to hash
	unless (open(REF, "$cancerGeneList")){
		print STDERR "Can not open file $cancerGeneList\n";
		exit;
	}
	my %CANCER_GENE;
	while(<REF>){
		chomp;
		$CANCER_GENE{$_} ='yes';
	}
	close REF;
	# HotSpot site 
	unless(open(HOT, "$hotspot")){
		print STDERR "Can not open file $hotspot\n";
		exit;
	}
	my %HOT_SPOT;
	while(<HOT>){
		chomp;
		my @local = split ("\t", $_);
		my $key = join "\t", @local[0..2];
		$HOT_SPOT{"$key"} = $local[3];
	}
	close HOT;
	# Germline mutation to work on.
	unless (open (ORI,"$germline")){
		print STDERR "Can not open file $germline\n";
		exit;
	}
	my $head =`grep -m1 -P "^Chr\tStart\tEnd\tRef\tAlt\t" $annotation |sort |uniq`;
	chomp($head);
	print "$head\tSample\tCaller\tQUAL\tFS\tTotalReads\tAltReads\tVAF\tSource\tLevel\n";
	while (<ORI>){
		chomp;
		my @temp = split("\t", $_);
		my $vcf;
		my $end = @temp - 1 ;
		my $site_sample = "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]";
		my $site = "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]";
		$vcf = join "\t", @temp[5..$end];
		if (!exists $DATA{$site}){ # i.e. position is germline!!
		#if (!exists $DATA{$site_sample}){ # i.e. position in sample is germline
			my ($source, $level) = findSource($ANNOTATION{"$site"});
			my @ANN = split("\t", $ANNOTATION{"$site"});
			my $vaf = VAF($temp[9], $temp[10]);
			if(($temp[9] !~ /\D/) and $temp[9] >=10 and $vaf >=0.25){
				if($source =~ /[ACMG-clinvar|hgmd|ACMG]/){
					if(exists $CANCER_GENE{$ANN[1]}){
						if ($ANN[0] =~ /splicing/ or $ANN[3] =~ /stopgain/ or $ANN[3]=~ /^frameshift/){
							$level = "stringent";
						}
						$source = $source.";CancerGene";
						if (exists $HOT_SPOT{"$temp[0]\t$temp[1]\t$temp[2]"}){
							$level  = "stringent";
							$source = $source.";".$HOT_SPOT{"$temp[0]\t$temp[1]\t$temp[2]"};
						}
					}
					else{
						if (exists $HOT_SPOT{"$temp[0]\t$temp[1]\t$temp[2]"}){
							$level  = "stringent";
							$source = $source.";".$HOT_SPOT{"$temp[0]\t$temp[1]\t$temp[2]"};
						}	
					}
					print "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$ANNOTATION{$site}\t$vcf\t$vaf\t$source\t$level\n";
				}
				elsif(exists $CANCER_GENE{$ANN[1]}){
					$level="2";
					if ($ANN[0] =~ /splicing/ or $ANN[3] =~ /stopgain/ or $ANN[3]=~ /^frameshift/){
						$level = "stringent";
					}
					$source = $source."CancerGene";
					print "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$ANNOTATION{$site}\t$vcf\t$vaf\t$source\t$level\n";
				}
			}
		}
	}
	close ORI;
}
sub findSource{
	my ($input)= (@_);
	my %source;
	my $level;
	my @ANN = split("\t", $input);
	if ($ANN[1] eq $ANN[147]){
		$source{'ACMG'} = 'yes';
		$level="2";
		if ($ANN[57] =~ /^Pathogenic/ or $ANN[57] =~ /\|Pathogenic/ or $ANN[57] =~ /^Likely Pathogenic/ or $ANN[57] =~ /\|Likely Pathogenic/){
			$source{'ACMG-clinvar'} = 'yes';
			$level = "stringent";
		}
	}
	if($ANN[64] =~ /^Disease causing mutation$/){  # HGMD
		$source{'HGMD'} = 'yes';
		$level="2";
		if ($ANN[57] =~ /^Pathogenic/ or $ANN[57] =~ /\|Pathogenic/ or $ANN[57] =~ /^Likely Pathogenic/ or $ANN[57] =~ /\|Likely Pathogenic/){
			$source{'HGMD-clinvar'} = 'yes';
			$level = "stringent";
		}
	}
	my $return = join(";", (sort keys %source));
	return($return, $level);
}
sub Somatic{
#/data/Clinomics/Ref/annovar/hg19_SomaticActionableSites.txt NCI0276/NCI0276/db/NCI0276.somatic
	my ($ref, $cgc ,$subject, $annotation) = (@_);
	unless (open(ANN_FH, "$ref")){
		print STDERR "\n\nCan not open $ref\n"; 
		exit;
	}
	my %DATA;
	while(<ANN_FH>){
		chomp;
		my @local = split ("\t", $_);
		my $key = join "\t", @local[0..2];
		my $end = @local - 1 ;
		my $value = join "\t", @local[3..$end];
		$DATA{"$key"} = $value;
	}
	close ANN_FH;
	unless (open(ANN_FH, "$cgc")){
		print STDERR "\n\nCan not open $ref\n";
		exit;
	}
	my %CGC;
	while(<ANN_FH>){
		chomp;
		my @local = split ("\t", $_);
		$CGC{"$local[0]"} = "CancerGeneCensus";
	}
	close ANN_FH;
	unless (open (ORI,"$subject")){
		print STDERR "\n\nCan not open $subject\n"; 
		exit;
	}
	#Annotations to hash
	unless (open(REF, "$annotation")){
		print STDERR "Can not open file $annotation\n";
		exit;
	}
	my %ANNOTATION;
	while(<REF>){
		chomp;
		my @local = split ("\t", $_);
		my $key = join "\t", @local[0..4];
		my $end = @local - 1 ;
		my $value = join "\t", @local[5..$end];
		$ANNOTATION{$key} =$value;
	}
	close REF;
	print "Chr\tStart\tEnd\tRef\tAlt\t";
	print $ANNOTATION{"Chr\tStart\tEnd\tRef\tAlt"};
	print "\tSample\tCaller\tQUAL\tFS\tTotalReads\tAltReads\tVAF\tSource\n";
	while (<ORI>){
		chomp;
		my @temp = split("\t", $_);
		my $val;
		my $vcf;	
		my $end = @temp - 1 ;
		my $vaf = VAF($temp[9], $temp[10]);
		my $key = join "\t", @temp[0..4];
		$val = "$temp[0]\t$temp[1]\t$temp[2]";
		$vcf = join "\t", @temp[5..$end];
		if (exists $DATA{$val}){
			print "$key\t$ANNOTATION{$key}\t$vcf\t$vaf\t$DATA{$val}\n";
		}
		else{
			my @local = split("\t",$ANNOTATION{$key});
			if (exists $CGC{$local[1]} and( $local[3] =~ /stopgain/ or $local[3]=~ /^frameshift/ or $local[0] =~ /splicing/)){
				print "$key\t$ANNOTATION{$key}\t$vcf\t$vaf\t$CGC{$local[1]}\n";
			}
		}
	}
	close ORI;
}

sub VAF{
	my ($total, $var) = (@_);
	my $vaf =0;
	if($var =~ /,/ or $total =~ /\./ or $total =~ /NA/){ 
		return($vaf);
	}
	elsif($total == 0){
		return($vaf);
	}
	else{
		$vaf = sprintf("%0.2f", $var/$total);
		return($vaf);
	}

}
