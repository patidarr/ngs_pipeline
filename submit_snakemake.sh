#!/bin/sh
#PBS -l walltime=96:00:00
#PBS -N ngs-pipeline
#
# Author: Rajesh Patidar
# 
# Usually slurm can translate the PBS varibles so no need to initialize the following sbatch vars.
#
##SBATCH --job-name="KhanLab"
##SBATCH --mail-type=FAIL
##SBATCH --cpus-per-task=1
##SBATCH --mem=1g
##SBATCH --gres=lscratch:01
#
TIME=$(date +"%Y_%m_%d")
#TIME=$(date +"%Y%m%d_%H")
if [[ `hostname` =~ "cn" ]] || [ `hostname` == 'biowulf.nih.gov' ]; then
	module use /data/khanlab/apps/Modules
	module load python/3.4.3
	export NGS_PIPELINE="/gpfs/gsfs4/users/khanlab/projects/patidar/ngs_pipeline/"
	export WORK_DIR="/gpfs/gsfs4/users/khanlab/projects/patidar/ngs_pipeline/test"
	export DATA_DIR="/data/khanlab/projects/DATA/"
	export ACT_DIR="/Actionable/"
	export HOST="biowulf.nih.gov"
	SNAKEFILE=$NGS_PIPELINE/ngs_pipeline.snakefile
	SAM_CONFIG=/gpfs/gsfs4/users/khanlab/projects/patidar/ngs_pipeline/samplesheet.json
elif [[ `hostname` =~ "tgcompute" ]] || [ `hostname` == 'login01' ] ; then
	module load python/3.4.3
	module load snakemake
	export NGS_PIPELINE="/projects/Clinomics/Tools/ngs_pipeline-dev"
	export WORK_DIR="/projects/Clinomics/Test_Run2/"
	export DATA_DIR="/projects/Clinomics/DATA/"
	export ACT_DIR="/Actionable/"
	export HOST="login01"
	SNAKEFILE=$NGS_PIPELINE/ngs_pipeline.snakefile
	SAM_CONFIG=$WORK_DIR/samplesheet.json
else 
	echo -e "Host `hostname` is not recognized\n"
	echo -e "This pipeline is customized to run on biowulf.nih.gov or TGen Cluster@ KhanLab\n";
	exit;
fi

cd $WORK_DIR
#if [ `cat $SAM_CONFIG |/usr/bin/json_verify -c` -ne "JSON is valid" ]; then
#       echo "$SAM_CONFIG is not a valid json file"
#       exit
#fi
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
	snakemake $cmd --cluster "sbatch --mail-type=FAIL -o log/{params.rulename}.%j.o {params.batch}" >& ngs_pipeline_${TIME}.log
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
