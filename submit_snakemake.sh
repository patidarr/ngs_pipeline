#!/bin/sh
#PBS -l walltime=96:00:00
#PBS -N ngs-pipeline
#
# Author: Rajesh Patidar
# 
# Usually slurm can translate the PBS varibles so no need to initialize the sbatch variables.
set -eo pipefail
if [[ $time == 'd' ]]; then
	export TIME="20160415"
elif [[ $time == 'p' ]]; then
	export TIME=$(date +"%Y%m%d")
else
	echo -e "Can not run without knowing which mode you would like to set time up\n";
	exit;
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

NOW=$(date +"%Y%m%d%H")
#export TIME=$(date +"%Y%m%d%H")
if [[ `hostname` =~ "cn" ]] || [ `hostname` == 'biowulf.nih.gov' ]; then
	module load snakemake/3.5.5.n1
	export HOST="biowulf.nih.gov"
elif [[ `hostname` =~ "tghighmem" ]] || [[ `hostname` =~ "tgcompute" ]] || [ `hostname` == 'login01' ] ; then
	module load snakemake/3.5.5
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
#if [ ! -d annovar ]; then
#	mkdir annovar
#	touch annovar/AnnotationInput.final.txt
#fi

export ACT_DIR="/Actionable/"
SNAKEFILE=$NGS_PIPELINE/ngs_pipeline.snakefile
SAM_CONFIG=$WORK_DIR/samplesheet.json

cmd="--directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --jobname {params.rulename}.{jobid} --nolock  --ri -k -p -T -j 3000 --stats ngs_pipeline_${NOW}.stats"
if [ $HOST   == 'biowulf.nih.gov' ]; then
	echo "Host identified as $HOST"
	snakemake $cmd --cluster "sbatch -o log/{params.rulename}.%j.o -e log/{params.rulename}.%j.e {params.batch}" >& ngs_pipeline_${NOW}.log
elif [ $HOST == 'login01' ]; then
	 echo "Host identified as $HOST"
	snakemake $cmd --cluster "qsub -W umask=022 -V -e $WORK_DIR/log/ -o $WORK_DIR/log/ {params.batch}" >& ngs_pipeline_${NOW}.log
fi
chmod -R 777 $WORK_DIR/.snakemake
