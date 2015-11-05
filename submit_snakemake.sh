#!/bin/sh
#SBATCH --job-name="Exome-Pipeline"
#SBATCH --mail-type=FAIL
#SBATCH --time="10-00:00:00"
#SBATCH --output=log/snakemake.%j.o
#SBATCH --partition="ccr"


#NOW=$(date +"%H%M%S_%m%d%Y")
NOW=$(date +"%Y%m%d_%H")

module load python/3.4.3
export NGS_PIPELINE="/data/khanlab/projects/patidar/Snakemake"
export WORK_DIR="/data/khanlab/projects/patidar/Snakemake"
SNAKEFILE=$NGS_PIPELINE/ngs_pipeline.rules
SAM_CONFIG=$WORK_DIR/samplesheet.json

if [ `cat $SAM_CONFIG |/usr/bin/json_verify -c` -ne "JSON is valid" ]; then
	echo "$SAM_CONFIG is not a valid json file"
	exit
fi

cd $WORK_DIR
if [ ! -d log ]; then
	mkdir log
fi

snakemake\
	--directory $WORK_DIR \
	--snakefile $SNAKEFILE \
	--configfile $SAM_CONFIG \
	--jobname '{rulename}.{jobid}' \
	--nolock \
	-k -p -T \
	-j 3000 \
	--stats ngs_pipeline_${NOW}.stats \
	--cluster "sbatch -o log/{params.rulename}.%j.o {params.batch}"\
	>& serpentine_${NOW}.log

# Summary 
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --summary

## DRY Run with Print out the shell commands that will be executed
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --dryrun -p

#DAG 
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --dag | dot -Tpng > dag.png

#Rulegraph
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG -n --forceall --rulegraph | dot -Tpng > rulegraph.png

# Mail Rulegraph and DAG to self
#  echo DAG |mutt -s "DAG" -a dag.png -a rulegraph.png -- patidarr@mail.nih.gov
