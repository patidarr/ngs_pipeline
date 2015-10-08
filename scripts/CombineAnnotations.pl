#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
###################################
# Author: Rajesh Patidar (rajbtpatidar@gmail.com)
# This takes a file list. and generate a Variant Matrix
# 	every file is in same variant format where last column have value for that file (article in this case)
# 
###################################

unless (open(LIST, $ARGV[0])){
	print STDERR "Please give a file containing list of pediatric studies;";
	exit;
}
chomp(my @files = <LIST>);
close LIST;
my %ALL_SITES;
my @files2 = @files;
foreach my $file (@files){ # Make the union of all the keys 
	print STDERR "Processing $file\n";
	fillHash($file); 
}

print STDERR "Created big hash\n";
# Add  column for every file. defaule space holder is "0";
foreach my $file (@files){
	my %TMP_HASH;
	print STDERR "Adding $file to big hash\n";
	unless (open (TMP, $file)){
		print STDERR "Can not open file $file\n";
		exit;
	}
	while(<TMP>){ # Read Files 
		my $line = $_;
		chomp $line;
		my @arr = split("\t", $line);
		my $key = "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]";
		$TMP_HASH{"$key"} = "$arr[5]";
	}
	close TMP;

	# Add new file to the big HASH
	foreach my $key (keys %ALL_SITES){
		if (exists $TMP_HASH{$key}){
			$ALL_SITES{$key} = "$ALL_SITES{$key}$TMP_HASH{$key},";
			delete $TMP_HASH{$key};
		}
		else{
			$ALL_SITES{$key} = "$ALL_SITES{$key}0,"; # Change default placeholder here.
		}
		
	}
	
}
# Final printing and some manipulations
foreach my $key (sort keys %ALL_SITES){
	chop($ALL_SITES{$key});
	my @values = split(",", $ALL_SITES{$key});
	if ($key =~ /^Chr/){
		print "$key\t".join("\t", @values)."\tPCG_TOTAL\tGRAND_TOTAL\n"; # 
		next;
	}
	my $pcg = sum(@values[1..$#values]);
	my $gt  = sum(@values);
	print "$key\t".join("\t", @values)."\t$pcg\t$gt\n";
}
###########################
### Make Big Hash
###########################
sub fillHash{
	my ($f) = @_;
	open(FH, $f);
	while(<FH>){
        	chomp;
        	my @temp = split ("\t", $_);
		my $key = "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]";
		$ALL_SITES{"$key"} = "";
	}
	close FH;
	return($f);
}
