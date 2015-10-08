#!/bin/sh
#BATCH --job-name="NCI0276"
#BATCH --time=10-00:00:00"
NOW=$(date +"%H%M%S_%m%d%Y")
module load python/3.4.3

WORK_DIR=/data/khanlab/projects/patidar/Testing/
SNAKEFILE=/data/khanlab/projects/patidar/Snakemake/ngs_pipeline.rules
SAM_CONFIG=/data/khanlab/projects/patidar/Snakemake/config_NCI0276.json

cd $WORK_DIR
snakemake\
	--directory $WORK_DIR \
	--snakefile $SNAKEFILE \
	--configfile $SAM_CONFIG \
	--jobname '{rulename}.{jobid}' \
	--nolock --notemp \
	-k -p -T \
	-j 3000 \
	--stats ngs_pipeline_${NOW}.stats \
	--cluster "sbatch -o log/{params.rulename}.%j.o {params.batch}" \
        >& ngs_pipeline_${NOW}.log





#--cluster "sbatch -e log/{params.rulename}.%j.e -o log/{params.rulename}.%j.o --partition={params.partition} --mem={params.mem} --time={params.time} {params.batch}"\
# Summary 
#snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --summary


## DRY Run with Print out the shell commands that will be executed
#snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --dryrun -p

#For saving this to a file
#snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --dag | dot -Tpdf > dag.pdf
#snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --dag | dot -Tpng > ~/dag.png

#snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG -n --forceall --rulegraph | dot -Tpng > ~/dag_rajesh.png
#snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG -n --forceall --rulegraph | dot -Tpdf > ~/dag_rajesh.pdf


#snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --forceall --dag | dot -Tpdf > dag.pdf
#snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --forceall --dag | dot -Tpng > dag.png

#echo DAG |mutt -s "DAG" -a dag.pdf -- patidarr@mail.nih.gov
