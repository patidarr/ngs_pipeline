#!/bin/sh
# This extract reads from defuse temp folder
# Author Rajesh Patidar
#sh extractReadsdeFuse.sh /data/khanlab/ref/GATK/hg19/defuse.config_hg19_ens69.txt temp raj defuse.filtered.txt
config=$1
outdir=$2
dumpdir=$3
defuseResult=$4

function nrwait() {
	local nrwait_my_arg
	if [[ -z $1 ]] ; then
		nrwait_my_arg=2
	else
		nrwait_my_arg=$1
	fi
	while [[ $(jobs -p | wc -l) -ge $nrwait_my_arg ]] ; 
	do
		sleep 0.33;
	done
}
for cluster in `cut -f 1 $defuseResult |grep -v cluster`; 
	do
		get_reads.pl -c $config -o $outdir -i $cluster >$dumpdir/$cluster.txt &
		nrwait $SLURM_CPUS_ON_NODE
	done
	wait
