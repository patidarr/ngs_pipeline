#!/usr/bin/perl
use strict;
use warnings;

my $FILE = $ARGV[0]; # Vcf like file Name
my $BAM = $ARGV[1];
unless (open(FH, "$FILE")){
	print STDERR "Can not find the file $FILE\n";
}
while(<FH>){
        chomp;
	my $line = $_;
        my @d = split("\t", $_);
	if($_ =~ /^Chr/ or $_ =~ /^#/){
		print "$line\tTotCov\tRefCov\tVarCov\tVAF\n";
	}
	else{
		my $rnaseq;		
		$d[0] =~ s/chr//g;
		$rnaseq = `samtools mpileup -d1000000000000 -r $d[0]:$d[1]-$d[1] "$BAM" 2>/scratch/out |cut -f 3-5`;
		chomp $rnaseq;
		if($rnaseq =~ /\d/){ # Have coverage 
			my @s = split("\t", $rnaseq);
			$s[2] =uc($s[2]);
			my $bases = join '', sort, split //, $s[2];
			my $Ref = () = ($bases =~ /$d[3]/g); # 4th Column is the Ref Column
			my $Alt = () = ($bases =~ /$d[4]/g); # 5th Column is the Alt Column 
			my $totalRef = $Ref + $Alt;
			if($totalRef >=1){ 
				if($Alt >=1){
					my $vaf = ($Alt/$totalRef);
					print "$line\t$totalRef\t$Ref\t$Alt\t$vaf\n";
				}
				else{
					print "$line\t$totalRef\t$Ref\t$Alt\t0\n";
				}
			}
			else{
				print "$line\t0\t0\t0\t0\n";
			}
		}
		else{ # no output from mpileup no COVERAGE
			print "$line\t0\t0\t0\t0\n";
		}
	}
}	
close FH;
