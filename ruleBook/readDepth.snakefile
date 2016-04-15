############
##      Depth Per Base for entire target intervals
############
rule ReadDepth:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.bwa.final.bam",
		bai="{base}/{TIME}/{sample}/{sample}.bwa.final.bam.bai"
	output:
		"{base}/{TIME}/{sample}/qc/{sample}.depth_per_base"
	version: config['bedtools']
	params:
		rulename	= "readDepth",
		target_intervals= lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]],
		batch		= config[config['host']]["job_bedtools"]
	shell: """
	#######################	
	module load bedtools/{version}
	module load R
	echo -e "chr\\tstart\\tend\\tgene\\tposition\\tdepth" >  {output}
	cut -f1-4 {params.target_intervals} | bedtools coverage -abam {input.bam} -b - -d >> {output}
	#######################
	"""
