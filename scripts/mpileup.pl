#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
no if $] >= 5.017011, warnings => 'experimental::smartmatch';
use 5.010;
local $SIG{__WARN__} = sub {
        my $message =shift;
        die $message;
};
my $FILE = $ARGV[0]; # Vcf like file Name
my $BAM  = $ARGV[1];
my $RNA  = $ARGV[2]; #RNASeq Result
unless (open(FH, "$FILE")){
	print STDERR "Can not find the file $FILE\n";
}
my @a = split(/[.]/, basename($BAM));
my %Expressed;
unless (open(IN2, "$RNA")){
        print STDERR "Can not find the file $FILE\n";
}
while(<IN2>){
	chomp;
	my @line = split("\t",$_);
	$Expressed{join("\t",@line[0..4])} = join("\t",@line[9..13]);
}
close IN2;
while(<FH>){
	chomp;
	my $line = $_;
	my @d = split("\t", $_);
	if($_ =~ /^Chr/ or $_ =~ /^#/){
		print "$line\t$a[0].GT\tRNASeq.TotCov\tRNASeq.RefCov\tRNASeq.VarCov\tRNASeq.VAF\n";
	}
	elsif ($d[1] eq $d[2] and length($d[3]) eq length($d[4]) and length($d[3]) <2){
		my $rnaseq;		
		$rnaseq = `samtools mpileup -d1000000000000 -r $d[0]:$d[1]-$d[1] "$BAM" 2>/dev/null |cut -f 3-5`;
		chomp $rnaseq;
		if($rnaseq =~ /\d/){ # Have coverage 
			my @s = split("\t", $rnaseq);
			if ($s[1] <1){
				print "$line\tNA\t0\t0\t0\t0\n";
			}
			else{
				$s[2] =uc($s[2]);
				my $bases = join '', sort, split //, $s[2];
				my $Ref = () = ($bases =~ /$d[3]/g); # 4th Column is the Ref Column
				my $Alt = () = ($bases =~ /$d[4]/g); # 5th Column is the Alt Column 
				my $totalRef = $Ref + $Alt;
				if($totalRef >=1){ 
					if($Alt >=1){
						my $vaf = ($Alt/$totalRef);
						print "$line\tNA\t$totalRef\t$Ref\t$Alt\t$vaf\n";
					}
					else{
						print "$line\tNA\t$totalRef\t$Ref\t$Alt\t0\n";
					}
				}
				else{
					print "$line\tNA\t0\t0\t0\t0\n";
				}
			}
		}
		else{ # no output from mpileup no COVERAGE
			print "$line\tNA\t0\t0\t0\t0\n";
		}
	}
	else{
		if (exists $Expressed{join("\t",@d[0..4])}){
			print "$_\t".$Expressed{join("\t",@d[0..4])}."\n";
		}
		else{
			my ($status,$info) = search($d[0], $d[1], $d[2]);
			if ($status eq 'match'){
				print "$_\t$info\n";
			}
			else{
				print "$_\tNA\t0\t0\t0\t0\n";

			}
		}

	}
}	
close FH;
sub search{
	my ($chr, $start, $end) = (@_);
	foreach my $indel (keys %Expressed){
		my @ind = split("\t", $indel);
		if ($ind[0] eq $chr and ($start ~~ [$ind[1] - 10 .. $ind[1] + 10]  or $end ~~ [$ind[2] - 10 .. $ind[2] + 10])){
			return("match", $Expressed{$indel});
		}
	}
}

