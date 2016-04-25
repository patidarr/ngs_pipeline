#!/bin/sh
#PBS -l walltime=96:00:00
#PBS -N ngs-pipeline
#
# Author: Rajesh Patidar
# 
# Usually slurm can translate the PBS varibles so no need to initialize the following sbatch vars.
set -eo pipefail
if [[ $time == 'd' ]]; then
	export TIME="20160415"
elif [[ $time == 'p' ]]; then
	export TIME=$(date +"%Y%m%d")
else
	echo -e "Can not run without knowing which mode you would like to set time up\n";
	exit;
fi
#export TIME=$(date +"%Y%m%d")
#export TIME=$(date +"%Y%m%d%H")
if [[ `hostname` =~ "cn" ]] || [ `hostname` == 'biowulf.nih.gov' ]; then
	module use /data/khanlab/apps/Modules
	module load python/3.4.3
	export NGS_PIPELINE="/data/khanlab/projects/patidar/ngs_pipeline/"
	export WORK_DIR="/data/khanlab/projects/DNASeq/"
	export DATA_DIR="/data/khanlab/projects/DATA/"
	export ACT_DIR="/Actionable/"
	export HOST="biowulf.nih.gov"
	SNAKEFILE=$NGS_PIPELINE/ngs_pipeline.snakefile
	SAM_CONFIG=$WORK_DIR/samplesheet.json
	
elif [[ `hostname` =~ "tghighmem" ]] || [[ `hostname` =~ "tgcompute" ]] || [ `hostname` == 'login01' ] ; then
	module use /home/patidarr/Modules
	module load python/3.4.3
	module load snakemake
	export NGS_PIPELINE="/projects/Clinomics/Tools/ngs_pipeline/"
	export WORK_DIR="/projects/Clinomics/Test_Run3/"
	export DATA_DIR="/projects/Clinomics/DATA/"
	export ACT_DIR="/Actionable/"
	export HOST="login01"
	SNAKEFILE=$NGS_PIPELINE/ngs_pipeline.snakefile
	SAM_CONFIG=$WORK_DIR/samplesheet.json
else 
	echo -e "Host `hostname` is not recognized\n"
	echo -e "This pipeline is customized to run on biowulf.nih.gov or TGen Cluster @ KhanLab\n";
	echo -e "If you would like to use it on another system, you have to change config/config_cluster.json and some hardcoded system dependent\n";
	exit;
fi

if [[ ! -z $dataDir ]]; then 
	export DATA_DIR=$dataDir
fi
if [[ ! -z $workDir ]]; then
	export WORK_DIR=$workDir
	SAM_CONFIG=$WORK_DIR/samplesheet.json
fi



cd $WORK_DIR
if [ ! -d log ]; then
	mkdir log
fi
if [ ! -d annovar ]; then
	mkdir annovar
	touch annovar/AnnotationInput.final.txt
fi
cmd="--directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --jobname {params.rulename}.{jobid} --nolock  -k -p -T -j 3000 --stats ngs_pipeline_${TIME}.stats"
if [ $HOST   == 'biowulf.nih.gov' ]; then
	echo "Host identified as $HOST"
	snakemake $cmd --cluster "sbatch -o log/{params.rulename}.%j.o -e log/{params.rulename}.%j.e {params.batch}" >& ngs_pipeline_${TIME}.log
elif [ $HOST == 'login01' ]; then
	 echo "Host identified as $HOST"
	snakemake $cmd --cluster "qsub  -V -e $WORK_DIR/log/ -o $WORK_DIR/log/ {params.batch}" >& ngs_pipeline_${TIME}.log
fi

# Summary 
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --summary

## DRY Run with Print out the shell commands that will be executed
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --dryrun -p -r

#DAG 
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --dag | dot -Tpng > dag.png

#Rulegraph
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG -n --forceall --rulegraph | dot -Tpng > Rulegraph.png

# Mail Rulegraph and DAG to self#  echo DAG |mutt -s "DAG" -a dag.png -a rulegraph.png -- patidarr@mail.nih.gov
#  echo DAG |mutt -s "DAG" -a dag.png -a Rulegraph.png -- patidarr@mail.nih.gov
