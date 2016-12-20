#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
local $SIG{__WARN__} = sub {my $message =shift; die $message;};
use Getopt::Long;
use Pod::Usage;
use List::Util qw(sum);
use File::Basename;
#################################
# This takes a comma seperated list of vcf files and make consensus vcf file.
#
#################################
my $VCF = '';
my $order='';
my $skip='';
my $help;
GetOptions(
		'vcf=s'		=>\$VCF,
		'order=s'	=>\$order,
		'filter=s'	=>\$skip,
		'help|h'	=>\$help,
	  )or pod2usage();

$help and pod2usage ();


if (!$VCF){
	print STDERR "Comma seperated list of vcf files required\n";
	exit;
}
my @listVCF=split(",",$VCF);
if ($#listVCF <1){
	print STDERR "You need 2 or more files\n";
	exit;
}
if (!$order){
        print STDERR "Comma seperated order of sample is required\n";
        exit;
}
my @samples=split(",", $order);
PrintHeader($listVCF[1]);
my @header;
sub PrintHeader{
	my ($file)=@_;
	unless (open(FH, $file)){
		print STDERR "\n\nCan not find file $file\n\n";
		exit;
	}
	while(<FH>){
		chomp;
		if ($_ =~ /^##/){
			print "$_\n";
		}
		if($_ =~ /^#CHROM/){
			my @h = split("\t", $_);
			my @i = @h[9..$#h];
			if(scalar(@i) ne scalar(@samples)){
				print STDERR "There are ".scalar(@i)." samples in vcf file but you provided ".scalar(@samples)." to order.\n";
				exit;
			}
			else{
				push @header, @samples ;
				print join("\t", @h[0..8],@samples)."\n";
			}
		}
	}
}
my %HASH;
#$idx_A = first { $format[$_] eq 'AU' } 0..$#format;
foreach my $vcf(@listVCF){
	unless (open(FH, $vcf)){
		print STDERR "\n\nCan not find file $vcf\n\n";
		exit;
	}
	my %ord;
#	print "Processing... $vcf\n";
	while(<FH>){
		chomp;
		if ($_ =~ /^##/){
			next;
		}
		if ($_ =~ /^#CHROM/){
			my @h = split("\t", $_);
			for (my $i=9; $i <= $#h; $i++){
				$ord{$h[$i]}  = $i;	
			}
			next;
		}

		my @line = split("\t", $_);
		if ($line[6] =~ /^$skip$/){
			next;
		}
		my $string = $line[8];
		foreach my $aa(@samples){
			$string .= "\t$line[$ord{$aa}]";
		}
		$HASH{"$line[0]\t$line[1]\t.\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t."} = "$string";
		
	}
}

foreach (sort keys %HASH){
	print "$_\t$HASH{$_}\n";
}

=head1 SYNOPSIS



 $0 -vcf 4060FrTu.MuTect.raw.vcf,4060FrTu.strelka.snvs.raw.vcf,4060FrTu.strelka.indels.raw.vcf -order 4060N,4060FrTu -filter REJECT
 This tool is to make merge vcf files generated on a set of samples using different tools e.g. somatic files generated on a pair of samples
	using MuTect,Strelka,VarScan, ect. Sometimes tools order the sampels in vcf files in either random order or in the order of bam files
	location on the command line, hence while merging the files order parameter is required.
 This keeps the record from the last file a record is seen. so if you like to keep the records from a particular file keep it at the end in vcf parameter
 filter parameter is a literal matching removal of lines from FILTER column of any vcf file the literal is found in. e.g. REJECT from MuTect vcf files.


 Usage:
	-h, -help, --help Print this message.
        -vcf    Comma seperated list of vcf files. (Required)
	-order	Comma seperated list of samples in vcf files (Required)
	-filter If some variants to be filtered out. Literal matching only. 

	For questions or comments, please contact: Rajesh Patidar <rajbtpatidar@gmail.com>
 Example usage:
	module load vcftools
	./consensusSomaticVCF.pl -vcf /data/khanlab/projects/processed_DATA/4060/20161207/4060FrTu/calls/4060FrTu.strelka.indels.raw.vcf,/data/khanlab/projects/processed_DATA/4060/20161207/4060FrTu/calls/4060FrTu.strelka.snvs.raw.vcf,/data/khanlab/projects/processed_DATA/4060/20161207/4060FrTu/calls/4060FrTu.MuTect.raw.vcf -order 4060N,4060FrTu -filter REJECT |vcf-subset -u -c 4060FrTu >4060FrTu.somatic.vcf
=cut
