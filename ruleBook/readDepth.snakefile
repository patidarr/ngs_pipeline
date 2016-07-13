############
##      Depth Per Base for entire target intervals
############
rule ReadDepth:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.bwa.final.bam",
		bai="{base}/{TIME}/{sample}/{sample}.bwa.final.bam.bai",
		target_intervals= lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]],
	output:
		"{base}/{TIME}/{sample}/qc/{sample}.depth_per_base"
	version: config['bedtools']
	params:
		rulename	= "readDepth",
		batch		= config[config['host']]["job_bedtools"]
	shell: """
	#######################	
	module load bedtools/{version}
	module load R
	echo -e "chr\\tstart\\tend\\tgene\\tposition\\tdepth" >  {output}
	cut -f1-4 {input.target_intervals} | bedtools coverage -abam {input.bam} -b - -d >> {output}
	#######################
	"""
