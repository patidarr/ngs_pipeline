############
# CoveragePlot
############
rule Ideogram:
	input:
		cnFile ="{subject}/{TIME}/{Tumor}/copyNumber/{Tumor}.copyNumber.txt",
		plot   =NGS_PIPELINE + "/scripts/ideogram.R",
		moving =NGS_PIPELINE + "/scripts/MovingWindow.pl", 
	output: 
		v1="{subject}/{TIME}/{Tumor}/copyNumber/{Tumor}.CN.raw.png",
		v2="{subject}/{TIME}/{Tumor}/copyNumber/{Tumor}.CN.png",
	version: config["R"]
	params:
		rulename = "ideogram",
		batch    = config[config['host']]["job_default"],
		genome   = config["reference"].replace('.fasta', '.genome'),
	shell: """
	#######################
	module load R
	module load bedtools
	grep -v Failed {input.cnFile} |awk '{{OFS="\t"}};{{print $1,$2,$3,$7,"-"}}' >${{LOCAL}}/{wildcards.sample}
	{input.plot} -f ${{LOCAL}}/{wildcards.sample} -o {output.v1} -t "{wildcards.Tumor}"

	bedtools makewindows -g {params.genome} -w 100000 | awk '{{OFS="\t"}}{{print $1,$2,$3,"win","end"}}' >${{LOCAL}}/{wildcards.sample}.movingWin.bed
	cat ${{LOCAL}}/{wildcards.sample} ${{LOCAL}}/{wildcards.sample}.movingWin.bed |sortBed -i - >${{LOCAL}}/{wildcards.sample}.sorted.bed
	perl {input.moving} ${{LOCAL}}/{wildcards.sample}.sorted.bed >${{LOCAL}}/{wildcards.sample}.sorted.moving.bed

	
	{input.plot} -f ${{LOCAL}}/{wildcards.sample}.sorted.moving.bed -o {output.v2} -t "{wildcards.Tumor}"
	cp {output.v1} {output.v2} {wildcards.subject}/{TIME}{ACT_DIR}/
	#######################
	"""
