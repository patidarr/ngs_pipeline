#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
no if $] >= 5.017011, warnings => 'experimental::smartmatch';
use 5.010;
local $SIG{__WARN__} = sub {my $message =shift; die $message;};
#tracking_id	class_code	nearest_ref_id	gene_id	gene_short_name	tss_id	locus	length	coverage	FPKM	FPKM_conf_lo	FPKM_conf_hi	FPKM_status
#TUBB8	-	-	TUBB8	TUBB8	-	10:92827-95178	-	-	0.457226	0.200609	0.713843	OK
#LOC101927762	-	-	LOC101927762	LOC101927762	-	10:978965-988683	-	-	0	0	0	OK
#DIP2C	-	-	DIP2C	DIP2C	-	10:320129-735608	-	-	0.707061	0.57992	0.834202	OK
#MIR7641-2	-	-	MIR7641-2	MIR7641-2	-	10:327976-328026	-	-	0	0	0	OK
#MIR5699	-	-	MIR5699	MIR5699	-	10:687628-687718	-	-	0	0	0	OK
#PRR26	-	-	PRR26	PRR26	-	10:695887-711109	-	-	0.0167573	0	0.0834863	OK
#LARP4B	-	-	LARP4B	LARP4B	-	10:852853-977645	-	-	6.96218	6.57729	7.34708	OK
#GTPBP4	-	-	GTPBP4	GTPBP4	-	10:1034348-1063708	-	-	40.4079	38.6503	42.1655	OK
#Tracking_ID	Gene_id	Gene_Symbol	tss_id	locus	log2(FPKM)
#LOC100506869	LOC100506869	LOC100506869	-	12:58985483-59206450	0.0000
#LOC102724604	LOC102724604	LOC102724604	-	3:181138950-181160274	0.0000
#MTVR2	MTVR2	MTVR2	-	17:54961462-54962274	0.0000

unless (open (FH, $ARGV[0])){
	print "can not open file $ARGV[0]\n";
	exit;
}
while(<FH>){
	chomp;
	my @line =split("\t", $_);
	if ($_=~ '^tracking_id'){
		print "Tracking_ID\tGene_ID\tGene_Name\tLocus\tlog2(FPKM)\n";
		next;
	}
	if ($line[6] =~ /^\d+:/ or $line[6] =~ /^[X|Y]:/){
		my $fpkm=$line[9] + 1;
		$fpkm=sprintf ("%.3f", log($fpkm)/log(2));
		print "$line[0]\t$line[3]\t$line[4]\t$line[6]\t$fpkm\n";
	}
}
close FH;
