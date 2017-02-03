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
		R 		= config['version_R'],
		samtools	= config['samtools'],
		batch		= config[config['host']]["job_bedtools"]
	shell: """
	#######################	
	module load bedtools/{version} samtools/{params.samtools}
	echo -e "chr\\tstart\\tend\\tgene\\tposition\\tdepth" >  {output}
	samtools view -hF 0x400 -q 30 -L {input.target_intervals} {input.bam} | samtools view -ShF 0x4 - | samtools view -SuF 0x200 - | bedtools coverage -abam - -b {input.target_intervals} -d >> {output}
	#cut -f1-4 {input.target_intervals} | bedtools coverage -abam {input.bam} -b - -d >> {output}
	#######################
	"""
