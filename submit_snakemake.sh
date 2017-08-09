#!/bin/sh
#PBS -N ngs-pipeline
#
# Author: Rajesh Patidar
# 
# Usually slurm can translate the PBS varibles so no need to initialize the sbatch variables.
#set -eo pipefail
module load snakemake/3.8.0
if [[ $time == 'd' ]]; then
	export TIME="20160415"
elif [[ $time == 'p' ]]; then
	export TIME=$(date +"%Y%m%d")
else
	export TIME=$time
#	echo -e "Can not run without knowing which mode you would like to set time up\n";
#	exit;
fi
if [[ ! -z $ngs ]]; then
        export NGS_PIPELINE=$ngs
fi
if [[ ! -z $dataDir ]]; then
        export DATA_DIR=$dataDir
fi

if [[ ! -z $workDir ]]; then
        export WORK_DIR=$workDir
        SAM_CONFIG=$WORK_DIR/samplesheet.json
fi

## for samplesheet
if [[ $sheet == 'samplesheet.json' ]]; then
	SAM_CONFIG=$WORK_DIR/samplesheet.json
else
	SAM_CONFIG=$sheet
fi
NOW=$(date +"%Y%m%d_%H%M%S")
#export TIME=$(date +"%Y%m%d%H")
export TMP="$NOW"
if [[ `hostname` =~ "cn" ]] || [ `hostname` == 'biowulf.nih.gov' ]; then
	#module load snakemake/3.5.5.nl
	export HOST="biowulf.nih.gov"
elif [[ `hostname` =~ "tghighmem" ]] || [[ `hostname` =~ "tgcompute" ]] || [ `hostname` == 'login01' ] ; then
	#module load snakemake/3.5.5
	export HOST="login01"
else 
	echo -e "Host `hostname` is not recognized\n"
	echo -e "This pipeline is customized to run on biowulf.nih.gov or TGen Cluster @ KhanLab\n";
	echo -e "If you would like to use it on another system, you have to change config/config_cluster.json and some hardcoded system dependencies\n";
	exit;
fi


cd $WORK_DIR
if [ ! -d log ]; then
	mkdir log
fi

export ACT_DIR="/Actionable/"
SNAKEFILE=$NGS_PIPELINE/ngs_pipeline.rules

cmd="--directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --jobscript $NGS_PIPELINE/scripts/jobscript.sh --jobname {params.rulename}.{jobid} --nolock  --ri -k -p -T -r -j 3000 --resources DeFuse=25 --resources SIFT=8 --stats ngs_pipeline_${NOW}.stats -R RNASeq "
#cmd="--directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --jobscript $NGS_PIPELINE/scripts/jobscript.sh --jobname {params.rulename}.{jobid} --nolock  --ri -k -p -T -r -j 3000 --resources DeFuse=25 --resources SIFT=8 --stats ngs_pipeline_${NOW}.stats -R makeConfig RNASeq "
umask 022
if [ $HOST   == 'biowulf.nih.gov' ]; then
	echo "Host identified as $HOST"
	echo "Variables are $cmd"
	snakemake $cmd --cluster "sbatch -o log/{params.rulename}.%j.o -e log/{params.rulename}.%j.e {params.batch}" >& ngs_pipeline_${NOW}.log
elif [ $HOST == 'login01' ]; then
	echo "Host identified as $HOST"
	echo "Variables are $cmd"
	snakemake $cmd --cluster "sbatchT -o log/{params.rulename}.%j.o -e log/{params.rulename}.%j.e {params.batch}" >& ngs_pipeline_${NOW}.log
	rm -rf /projects/scratch/ngs_pipeline_${NOW}_*
fi

if [ -f ngs_pipeline_${NOW}.stats ]; then
	python $NGS_PIPELINE/scripts/stats2Table.py ngs_pipeline_${NOW}.stats >ngs_pipeline_${NOW}.stats.txt
fi


if [ -s ngs_pipeline_${NOW}.stats.txt ]; then
	print "No job took significant time"
else
	rm -rf ngs_pipeline_${NOW}.stats.txt	
fi
