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
my $FILE = $ARGV[0]; # Hotspot input file
my $REF  = $ARGV[1];
my $BAM  = $ARGV[2]; # Bam file in question
my $sample=$ARGV[3];
my $type=$ARGV[4];
my $capture=$ARGV[5];
my $pileup;
$pileup = `samtools mpileup -d1000000000000 -Q 20 -q 30 -f $REF -l $FILE "$BAM" 2>/dev/null|cut -f 1,2,3,5`;
chomp $pileup;
my @input = split("\n", $pileup);
foreach my $new_line(@input){
	my @s = split("\t", $new_line);
	if (defined $s[3]){
		$s[3] =uc($s[3]);
		$s[3] =~ s/\$//g;
		$s[3] =~ s/\[//g;
		$s[3] =~ s/\]//g;
		$s[3] =~ s/\^//g;
		$s[3] =~ s/\./$s[2]/g;
		$s[3] =~ s/\,/$s[2]/g;
		my $bases = join '', sort, split //, $s[3];
		my $Ref = () = ($bases =~ /$s[2]/g);
		my $A = () = ($bases =~ /A/g);
		my $T = () = ($bases =~ /T/g);
		my $G = () = ($bases =~ /G/g);
		my $C = () = ($bases =~ /C/g);
		my $total = $A+$T+$G+$C;
		if ($total eq $Ref){
			next;
		}
		else{
			foreach my $base("A", "T","G","C"){
				if ($base eq $s[2]){
					next;
				}
				elsif((() = ($bases =~ /$base/g)) >3){
					my $alt = (() = ($bases =~ /$base/g));
					print "$s[0]\t$s[1]\t$s[1]\t$s[2]\t$base\t$sample\t$type\t$capture\tmpileup\tNA\tNA\t$total\t$alt\n";
				}
			}
		}
	}
}
