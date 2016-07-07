#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
use 5.010;
local $SIG{__WARN__} = sub {
        my $message =shift;
        die $message;
};
###################################
# Author: Rajesh Patidar (rajbtpatidar@gmail.com)
# Combine annotations 
###################################
unless (open(LIST, $ARGV[0])){
	print STDERR "Please give a file containing list of individual files studies;";
	exit;
}
chomp(my @files = <LIST>);
close LIST;
my %ALL_SITES;
unless (open(FH, "$files[0]")){
	print STDERR "$files[0]\n";
	die;
}
while(<FH>){
	chomp;
	my @temp = split ("\t", $_);
	my $key = "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]";
	$ALL_SITES{"$key"} = "";
}
close FH;
shift @files;
# Add  column for every file. defaule space holder is "0";
foreach my $file (@files){
	my $cols=0;
	my %TMP_HASH;
	unless (open (TMP, $file)){
		print STDERR "Can not open file $file\n";
		die;
	}
	while(<TMP>){ # Read Files 
		chomp $_;
		$_ =~ s/^#//g;
		my @arr = split("\t", $_);
		my $key = join "\t", @arr[0..4];
		my $end = @arr - 1 ;
		$cols = $end - 5 + 1;
		my $value = join "\t", @arr[5..$end];
		if (! exists$TMP_HASH{"$key"} ){
			$TMP_HASH{"$key"} = "$value";
		}
	}
	close TMP;
	my $string="0\t"x$cols;
	# Add new file to the big HASH
	foreach my $key (keys %ALL_SITES){
		if (exists $TMP_HASH{$key}){
			$ALL_SITES{$key} = "$ALL_SITES{$key}$TMP_HASH{$key}\t";
			delete $TMP_HASH{$key};
		}
		else{
			$ALL_SITES{$key} = "$ALL_SITES{$key}$string";
		}
	}
}
# Final printing and some manipulations
foreach my $key (sort keys %ALL_SITES){
	chop($ALL_SITES{$key});
	print "$key\t$ALL_SITES{$key}\n";
}
