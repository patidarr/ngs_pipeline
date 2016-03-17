#!/usr/bin/perl
use strict;
use warnings;
open(FH, $ARGV[0]);
my %NB;

my $total=`wc -l $ARGV[0]|sed -e 's/ /\t/g'|cut -f 1 `;
while(<FH>){
	chomp;
	my @d =split("\t");
	$NB{"$d[0]\t$d[1]"} = "$d[2]";
}
close FH;

open(FH, $ARGV[1]);
my $count=0;
my $count1=0;
while(<FH>){
	chomp;
	my @d =split("\t");
	if(exists $NB{"$d[0]\t$d[1]"}){
		$count++;
		if($d[2] eq $NB{"$d[0]\t$d[1]"}){
			$count1++;
		}
	}
}
close FH;
if ($count >0){
	$count=$count1 / $count;
	$count=sprintf("%.2f", $count);
	print "$count\n";
}
else{
	print "0\n";
}
