#!/bin/sh
#SBATCH --job-name="KhanLab"
#SBATCH --mail-type=FAIL
#SBATCH --output=log/snakemake.%j.o
#SBATCH --partition="unlimited"
#SBATCH --cpus-per-task=1
#SBATCH --mem=1g

#
# Make sure log directory exists or change the output file location on line 4.
#
#
#NOW=$(date +"%H%M%S_%m%d%Y")
NOW=$(date +"%Y%m%d_%H")
module use /data/khanlab/apps/Modules
module load python/3.4.3
export NGS_PIPELINE="/data/khanlab/projects/patidar/Snakemake"
export WORK_DIR="/data/khanlab/projects/patidar/Snakemake"
export DATA_DIR="/data/khanlab/projects/DATA"
SNAKEFILE=$NGS_PIPELINE/ngs_pipeline.rules
SAM_CONFIG=$WORK_DIR/samplesheet.json

if [ `cat $SAM_CONFIG |/usr/bin/json_verify -c` -ne "JSON is valid" ]; then
	echo "$SAM_CONFIG is not a valid json file"
	exit
fi

cd $WORK_DIR

snakemake\
	--directory $WORK_DIR \
	--snakefile $SNAKEFILE \
	--configfile $SAM_CONFIG \
	--jobname '{rulename}.{jobid}' \
	--nolock \
	-k -p -T \
	-j 3000 \
	--stats ngs_pipeline_${NOW}.stats \
	--cluster "sbatch --mail-type=FAIL -o log/{params.rulename}.%j.o {params.batch}"\
	>& ngs_pipeline_${NOW}.log

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
