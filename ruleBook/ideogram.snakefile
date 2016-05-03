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
	grep -v Failed {input.cnFile} |awk '{{OFS="\\t"}};{{print $1,$2,$3,$7,"-"}}' >${{LOCAL}}/{wildcards.Tumor}

	echo -e "CHR\\tPOS\\tP"  >${{LOCAL}}/{wildcards.Tumor}.temp
	grep -v Failed {input.cnFile} |grep -v "Chr" |sed -e 's/chr//g'|sed -e 's/X/23/g'|sed -e 's/Y/24/g' |sort -n|awk '{{OFS="\\t"}};{{print $1,($2+$3/2),$7}}' >>${{LOCAL}}/{wildcards.Tumor}.temp
	
	{input.plot} -f ${{LOCAL}}/{wildcards.Tumor}.temp -o {output.v1} -t "{wildcards.Tumor}"
	
	bedtools makewindows -g {params.genome} -w 100000 | awk '{{OFS="\\t"}}{{print $1,$2,$3,"win","end"}}' >${{LOCAL}}/{wildcards.Tumor}.movingWin.bed
	cat ${{LOCAL}}/{wildcards.Tumor} ${{LOCAL}}/{wildcards.Tumor}.movingWin.bed |sortBed -i - >${{LOCAL}}/{wildcards.Tumor}.sorted.bed
	perl {input.moving} ${{LOCAL}}/{wildcards.Tumor}.sorted.bed >${{LOCAL}}/{wildcards.Tumor}.sorted.moving.bed


	echo -e "CHR\\tPOS\\tP"  >${{LOCAL}}/{wildcards.Tumor}.sorted.moving.bed.temp
	cat ${{LOCAL}}/{wildcards.Tumor}.sorted.moving.bed |sed -e 's/chr//g'|sed -e 's/X/23/g'|sed -e 's/Y/24/g' |sort -n|awk '{{OFS="\\t"}};{{print $1,($2+$3/2),$4}}'>>${{LOCAL}}/{wildcards.Tumor}.sorted.moving.bed.temp
	{input.plot} -f ${{LOCAL}}/{wildcards.Tumor}.sorted.moving.bed.temp -o {output.v2} -t "{wildcards.Tumor}"
	cp {output.v1} {output.v2} {wildcards.subject}/{TIME}{ACT_DIR}/
	#######################
	"""
