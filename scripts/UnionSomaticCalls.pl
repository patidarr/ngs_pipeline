#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use List::Util qw(first);
use 5.010;
local $SIG{__WARN__} = sub {
	my $message =shift;
	die $message;
};
#
# Author Rajesh Patidar rajbtpatidar@gmail.com
#
#
# Strelka Expressed or Strelka (SNPS)
# Strelka Expressed or Strelka (Indels)
# muTect expressed or MuTect

my %CALLS;

foreach my $file(@ARGV){
	unless (open(FH, $file)){
		print STDERR "Can not open the file $file\n";
		exit;
	}
	my @a = split(/[.]/, basename($file));
	while(<FH>){
		chomp;
		my @line = split("\t", $_);
		my $key  = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]";
		if ($_ =~ /^Chr/){
			$CALLS{$key} = "$_\tCaller";
			next
		}
		$CALLS{$key} = "$_\t$a[1]";
	}
	close FH;
}
foreach (sort keys %CALLS) {
    print "$CALLS{$_}\n";
}
