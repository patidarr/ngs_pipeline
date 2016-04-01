#!/bin/sh
#SBATCH --job-name="KhanLab"
#SBATCH --mail-type=FAIL
#SBATCH --partition="unlimited"
#SBATCH --cpus-per-task=1
#SBATCH --mem=1g
# Make sure log directory exists or change the output file location on line 4.
#
#
#NOW=$(date +"%H%M%S_%m%d%Y")
NOW=$(date +"%Y%m%d_%H")
module use /data/khanlab/apps/Modules
module load python/3.4.3
module load snakemake
if [[ `hostname` =~ "cn" ]] || [ `hostname` == 'biowulf.nih.gov' ]; then
	export NGS_PIPELINE="/data/khanlab/projects/patidar/Snakemake"
	export WORK_DIR="/data/khanlab/projects/patidar/Snakemake/test"
	export DATA_DIR="/data/khanlab/projects/DATA"
	export ACT_DIR="/Actionable/"
	export HOST="biowulf.nih.gov"

elif [[ `hostname` =~ "tgcompute" ]] || [ `hostname` == 'login01' ] ; then
	export NGS_PIPELINE="/home/patidarr/ngs_pipeline/"
	export WORK_DIR="/home/patidarr/trial1/"
	export DATA_DIR="/projects/Clinomics/DATA/"
	export ACT_DIR="/Actionable/"
	export HOST="login01"
else 
	echo "Host `hostname` is not recognized\n"
	exit;
fi

SNAKEFILE=$NGS_PIPELINE/ngs_pipeline.rules
SAM_CONFIG=$WORK_DIR/samplesheet.json

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



cmd="--directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --jobname {params.rulename}.{jobid} --nolock  -k -p -T -j 3000 --stats ngs_pipeline_${NOW}.stats"
if [ $HOST   == 'biowulf.nih.gov' ]; then
	echo "Host identified as $HOST"
	snakemake $cmd --cluster "sbatch --mail-type=FAIL -o log/{params.rulename}.%j.o {params.batch}" >& ngs_pipeline_${NOW}.log
elif [ $HOST == 'login01' ]; then
	 echo "Host identified as $HOST"
	snakemake $cmd --cluster "qsub  -V -o $WORK_DIR/log/ {params.batch}" >& ngs_pipeline_${NOW}.log
fi


# Summary 
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --summary

## DRY Run with Print out the shell commands that will be executed
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --dryrun -p -r

#DAG 
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --dag | dot -Tpng > dag.png

#Rulegraph
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG -n --forceall --rulegraph | dot -Tpng > rulegraph.png

# Mail Rulegraph and DAG to self
#  echo DAG |mutt -s "DAG" -a dag.png -a rulegraph.png -- patidarr@mail.nih.gov
