#!/usr/bin/perl
use strict;
use warnings;
if ($ARGV[0] =~ 'filter'){
	filter();
}
elsif($ARGV[0] =~ 'collect'){
	collect();
}
else{
	print STDERR "$ARGV[0] does not looks like a valid argument\n";
	die;
}
if (!defined $ARGV[2]){
	print STDERR "Minimum is not defined $ARGV[2]\n";
	die;
}
if (!defined $ARGV[3]){
        print STDERR "Maximum is not defined $ARGV[3]\n";
        die;
}

sub filter{
	my $file= $ARGV[1];
	my $min=$ARGV[2];
	my $max=$ARGV[3];
	my $cgc=$ARGV[4];
	my $idx_col=6;
	my $idx_normal=3;
	open(FH, "$cgc");
        my %HASH;
        while(<FH>){
                chomp;
                my @a =split("\t", $_);
                $HASH{$a[0]} = $a[1],
        }
        close FH;
	
	open(FH, $file);
	while(<FH>){
		chomp;
		my @temp =split("\t", $_);
		my $status="normal";
		if ($temp[$idx_col] =~ /LRR/){
			print "$_\tStatus\tMatched Gene\tSource\n";
			next;
		}
		if($temp[$idx_normal] <30){
                        $status = "Failed";
                }
		elsif ($temp[$idx_col] <$min){
			$status = "Loss";
		}
		elsif($temp[$idx_col] >$max){
			$status = "Gain";	
		}
		my @GENES = split("; ", $temp[7]);
		foreach my $gene(@GENES){
			if (exists $HASH{$gene}){
				print "$_\t$status\t$gene\t$HASH{$gene}\n";
			}
			else{
				print "$_\t$status\t-\t-\n";
			}
		}
	}
}
sub collect{
	my $file =$ARGV[1];
	my $cgc  =$ARGV[2]; 
	my $hot  =$ARGV[3];
	my $all  =$ARGV[4];
	open(FH, "$cgc");
	my %HASH;
	while(<FH>){
		chomp;
		my @a =split("\t", $_);
		$HASH{$a[0]} = 'CancerGenesCensus',
	}
	close FH;
	open(FH, "$hot");
	while(<FH>){
		chomp;
		if (exists $HASH{$_}){
			$HASH{$_} = $HASH{$_}.', HotSpotGene';
		}
		else{
			$HASH{$_} = 'HotSpotGene';
		}
	}
	close FH;
	print "#Chr\tStart\tEnd\tNormalCoverage\tTumorCoverage\tRatio\tLRR\tGene(s)\tStrand(s)\tMatched Gene\tSource\n";
	my @list = `cut -f8 $file|sed -e "s/; /\\n/g" |sort |uniq`;
	chomp @list;
	for my $gene(@list){
		if (exists $HASH{$gene}){
			my $main = `grep -P "\t$gene\$|\t$gene;|; $gene;|; $gene\$" $all`;
			my @lines = split("\n", $main);
			for my $exon(@lines){
				print "$exon\t$gene\t$HASH{$gene}\n";
			}
		}
	}
}
