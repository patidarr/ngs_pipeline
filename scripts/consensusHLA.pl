#!/usr/bin/perl -ws
use 5.010;
local $SIG{__WARN__} = sub {my $message =shift; die $message;};
my $seq2HLA="$ARGV[0]";
my $HLAminer="$ARGV[1]";
print STDERR "consensusHLA.pl running on the $HLAminer and $seq2HLA\n";
#4060/20161207/4060N/HLA/HLAminer/HLAminer_HPTASR.csv
#4060/20161207/4060N/HLA/seq2HLA/4060N-ClassI.HLAgenotype4digits
open(IN,$seq2HLA);
my %HASH;
while(<IN>){
	chomp;
	next if ($_ =~ /#/);
	#Locus  Allele 1        Confidence      Allele 2        Confidence
	#A       A*33:01'        4.482428e-05    A*33:01 NA
	#B       B*51:01'        0.005186588     B*14:03'        0.01333522
	#C       C*02:02 0.001444297     C*08:02'        0.6296239
	#
	#
	$_ =~ s/'//g;
	my @line=split("\t",$_);
	if ($line[2] !~ /NA/i){
		$line[2] = sprintf("%.10f",$line[2]);
		$HASH{"$line[1]"} = $line[2] if (not exists $HASH{"$line[1]"});
	}
	if ($line[4] !~ /NA/i){
                $line[4] = sprintf("%.10f",$line[4]);
		$HASH{"$line[3]"} = $line[4] if (not exists $HASH{"$line[3]"});
	}
}
close IN;
print "#Allele\tseq2HLA_Confidence\tHLAminer_Score\n";
$/="------------------------------------------------------------------------------------------";
open(IN,$HLAminer);
while(<IN>){
	chomp;
	if ($. eq '3'){
		my @line=split("\n",$_);
		foreach my $ln (@line){
			$ln =~ s/\s+//g;
			if ($ln =~ /,/){
				my @allele = split(",",$ln);
				$allele[0] =~ s/P//;
				$allele[0] =~ s/Q//;
				if (length($allele[0]) >7){
					next;
				}
				if(exists $HASH{$allele[0]}){
					print "HLA-$allele[0]\t$HASH{$allele[0]}\t$allele[1]\n";
					delete $HASH{$allele[0]};
				}
				else{
					print "HLA-$allele[0]\tNotCalled\t$allele[1]\n";
				}
			}
		}
	}
	
}
foreach (keys %HASH){
	if($_ =~ /na/i){next;}
	print "HLA-$_\t$HASH{$_}\tNotCalled\n";
}
