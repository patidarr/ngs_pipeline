#!/bin/sh
# This extract reads from defuse temp folder 
# Author Rajesh Patidar
#sh extractReadsdeFuse.sh /data/khanlab/ref/GATK/hg19/defuse.config_hg19_ens69.txt temp raj defuse.filtered.txt
echo ${SLURM_CPUS_ON_NODE}
config=$1
outdir=$2
dumpdir=$3
defuseResult=$4

count=0
for cluster in `cut -f 1 $defuseResult |grep -v cluster`;
	do
		get_reads.pl -c $config -o $outdir -i $cluster >$dumpdir/$cluster.txt &
		let count+=1
		[[ $((count%SLURM_CPUS_ON_NODE)) -eq 0 ]] && wait
	done
