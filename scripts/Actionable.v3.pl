#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(first);
my $index_of_clinvar=57; # This does not include chr/alt, Func_RefGene==0
my $index_of_ACMG=182; 
my $index_of_HGMD=64;
my $index_of_Gene=1;
my $idx_anno_region=0;
my $idx_anno_eff=3;
if($ARGV[0] eq 'somatic'){
	Somatic($ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4], $ARGV[5]);
	# 1 == HotspotFile [reference]
	# 2 == CancerGene Census [reference]
	# 3 == combinedGeneList [reference]
	# 4 == somaticFile
	# 5 == Annotations
}
elsif($ARGV[0] eq 'germline' or $ARGV[0] eq 'rnaseq'){
	Germline($ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4], $ARGV[5]);
	# 1  == somaticFile
	# 2  == germlineFile
	# 3  == Annotations
	# 4  == CancerGene List (Inherited Diseases,JW_germline,ClinOmics Tier2,Genetics_HumanRef,ClinomicsPanel,FoundationMed,CancerGeneCensus,)
	# 5  == HotspotFile [reference]
}
else{
	die $!;
}
sub OpneFH{
	my ($file) =(@_);
	my $FH;
	unless (open($FH, "$file")){
		print STDERR "Can not open file $file\n";
		exit;
	}
	return($FH);
}
sub FillHASH{
	my ($FH) =(@_);
	my %HASH;
	while(<$FH>){
		chomp;
		$HASH{$_} ='yes';
	}
	close $FH;
	return(%HASH);
}
sub Germline{
	my ($somatic, $germline, $annotation, $cancerGeneList, $hotspot) = (@_);
	my $SOM = OpneFH($somatic);
	my $ANN = OpneFH($annotation);
	my $GL  = OpneFH($cancerGeneList);
	my $HOT = OpneFH($hotspot);
	my %DATA;
	while(<$SOM>){
		chomp;
		my @local = split ("\t", $_);
		my $key = join "\t", @local[0..4];
		$DATA{"$key"} = "somatic";
	}
	close $SOM;
	my %ANNOTATION;
	while(<$ANN>){
		chomp;
		my @local = split ("\t", $_);
		my $key = join "\t", @local[0..4];
		my $end = @local - 1 ;
		my $value = join "\t", @local[5..$end];
		$ANNOTATION{$key} =$value;
	}
	close $ANN;
	# Source from the combined file
	my %SOURCE;
	while(<$GL>){
		chomp;
		my @local = split ("\t", $_);
		$SOURCE{"$local[0]"} = "$local[1]";
	}
	close $GL;
	# HotSpot site 
	my %HOT_SPOT;
	while(<$HOT>){
		chomp;
		my @local = split ("\t", $_);
		my $key = join "\t", @local[0..2];
		$HOT_SPOT{"$key"} = $local[3];
	}
	close $HOT;
	# Germline mutation to work on.
	my $ORI =OpneFH($germline);
	my $head =`grep -m1 -P "^Chr\tStart\tEnd\tRef\tAlt\t" $annotation |sort |uniq`;
	chomp($head);
	print "$head\tSample\tSampleType\tCaptureType\tCaller\tQUAL\tFS\tTotalReads\tAltReads\tVAF\tSource\tLevel\n";
	my %Germline;
	while (<$ORI>){
		chomp;
		my @temp = split("\t", $_);
		my $vcf;
		my $end = @temp - 1 ;
		my $site = "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]";
		$vcf = join "\t", @temp[5..$end];
		if ($temp[6] =~ /RNASeq/){
			$Germline{$site} = $temp[7];
		}
		elsif($temp[6] =~ /Normal/){
			$Germline{$site} = $temp[7];
		}
		# Position is called as germline 
		# Position is not somatic called!!
		# in tumor the capture is same as in normal
		if (exists $Germline{$site} and !exists $DATA{$site} and $Germline{$site} eq $temp[7]){
			my %level;
			$level{"5"} = 'yes';
			my $source = "";
			my @ANN = split("\t", $ANNOTATION{"$site"});
			my $vaf = VAF($temp[11], $temp[12]);
			if($ANN[$index_of_HGMD] =~ /^Disease causing mutation$/){  # HGMD
                		if ($ANN[$index_of_clinvar] =~ /^Pathogenic/i or $ANN[$index_of_clinvar] =~ /\|Pathogenic/i or $ANN[$index_of_clinvar] =~ /^Likely Pathogenic/i or $ANN[$index_of_clinvar] =~ /\|Likely Pathogenic/i){
					$source = 'HGMD-clinvar';
					$level{"1"}  = "yes";
				}
				else{
					$source = 'HGMD';
					$level{"4"} = "yes";
				}
			}
			if (exists $SOURCE{$ANN[$index_of_Gene]}){ # ACMG
				$source =$SOURCE{$ANN[$index_of_Gene]};
				if ($ANN[$index_of_Gene] eq $ANN[$index_of_ACMG]){
					if ($ANN[$index_of_clinvar] =~ /^Pathogenic/i or $ANN[$index_of_clinvar] =~ /\|Pathogenic/i or $ANN[$index_of_clinvar] =~ /^Likely Pathogenic/i or $ANN[$index_of_clinvar] =~ /\|Likely Pathogenic/i){
						$level{"1"} = "yes";
					}
					else{
						$level{"3"} = "yes";
					}
				}
				if (exists $HOT_SPOT{"$temp[0]\t$temp[1]\t$temp[2]"}){
					$level{"1"} = "yes";
					$source = $source.";".$HOT_SPOT{"$temp[0]\t$temp[1]\t$temp[2]"};
				}
				if ($source =~ /^InheritedDiseases$/){
					if ($ANN[$idx_anno_region] =~ /splicing/ or $ANN[$idx_anno_eff] =~ /stopgain/ or $ANN[$idx_anno_eff]=~ /^frameshift/){
						$level{"3"} = "yes";
					}
					elsif($vaf >=0.75){
						$level{"3"} = "yes";
					}
					else{
						$level{"4"} = "yes";
					}
				}
				else{
					if ($ANN[$index_of_clinvar] =~ /^Pathogenic/i or $ANN[$index_of_clinvar] =~ /\|Pathogenic/i or $ANN[$index_of_clinvar] =~ /^Likely Pathogenic/i or $ANN[$index_of_clinvar] =~ /\|Likely Pathogenic/i){
						$level{"1"} = "yes";
					}
					if ($ANN[$idx_anno_region] =~ /splicing/ or $ANN[$idx_anno_eff] =~ /stopgain/ or $ANN[$idx_anno_eff]=~ /^frameshift/){
						$level{"1"} = "yes";
					}
					else{
						$level{"3"} = "yes";
					}
				}
			}
			my @l = sort keys %level;
			if ($l[0] <=4){
				print "$site\t$ANNOTATION{$site}\t$vcf\t$vaf\t$source\tTier$l[0]\n";
			}
		}
	}
	close $ORI;
}
sub Somatic{
#/data/Clinomics/Ref/annovar/hg19_SomaticActionableSites.txt NCI0276/NCI0276/db/NCI0276.somatic
	my ($ref, $cgc, $geneList ,$subject, $annotation) = (@_);
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
	unless (open(ANN_FH, "$geneList")){
		print STDERR "\n\nCan not open $ref\n";
		exit;
	}
	my %CGC;
	while(<ANN_FH>){
		chomp;
		my @local  = split("\t", $_); 
		$CGC{"$local[0]"} = "$local[1]";
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
	print "\tSample\tSampleType\tCaptureType\tCaller\tQUAL\tFS\tTotalReads\tAltReads\tVAF\tSource\tLevel\n";
	while (<ORI>){
		chomp;
		my @temp = split("\t", $_);
		my $val;
		my $vcf;
		my $end = @temp - 1 ;
		my $vaf = VAF($temp[11], $temp[12]);
		my $key = join "\t", @temp[0..4];
		$val = "$temp[0]\t$temp[1]\t$temp[2]";
		$vcf = join "\t", @temp[5..$end];
		if (exists $DATA{$val}){
			print "$key\t$ANNOTATION{$key}\t$vcf\t$vaf\t$DATA{$val}\tTier1\n";
		}
		else{
			my @local = split("\t",$ANNOTATION{$key});
			if (exists $CGC{$local[1]}){
				if ($local[3] =~ /stopgain/ or $local[3]=~ /^frameshift/ or $local[0] =~ /splicing/){
					print "$key\t$ANNOTATION{$key}\t$vcf\t$vaf\t$CGC{$local[1]}\tTier2\n";
				}
				else{
					print "$key\t$ANNOTATION{$key}\t$vcf\t$vaf\t$CGC{$local[1]}\tTier3\n";
				}
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
