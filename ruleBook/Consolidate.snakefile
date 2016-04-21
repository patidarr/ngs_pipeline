rule QC:
	input: bam="{base}/{TIME}/{sample}/{sample}.bwa.final.bam"
	output: "{base}/{TIME}/{sample}/qc/{sample}.consolidated_QC"
	version:
		config['samtools']
	params:
		rulename	= "consolidated_QC",
		target_intervals= lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]],
		bedtools	= config['bedtools'],
		tool		= NGS_PIPELINE+ "/scripts/QC_stats_Final.py",
		batch   	= config[config['host']]["job_c_QC"],
		diagnosis	= config['Diagnosis'][sample]
	shell:	"""
	#######################
	module load python/2.7.9
	module load samtools/{version}
	module load bedtools/{params.bedtools}
	python {params.tool} {input.bam} {params.target_intervals} ${{LOCAL}} {wildcards.base} {wildcards.sample}  "{params.diagnosis}" > {output}
	#######################
	"""
rule QC_Summary:
	input : CON_QC
	output: "Consolidated_QC.txt"
	params:
		rulename  = "QC_Sum",
		batch	  = config[config['host']]["job_default"]
	shell: """
	#######################
	cat {input} |sort |uniq >{output}
	#######################
	"""
