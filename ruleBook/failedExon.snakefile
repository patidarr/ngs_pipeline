TARGET     += ["{subject}/{TIME}/{sample}/qc/{sample}.failExons".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
TARGET     += ["{subject}/{TIME}/{sample}/qc/{sample}.failGenes".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
rule  FailedExons:
	input: depth_file="{base}/{TIME}/{sample}/qc/{sample}.depth_per_base"
	output:
		failed_exon_stats="{base}/{TIME}/{sample}/qc/{sample}.failExons",
		failed_gene_stats="{base}/{TIME}/{sample}/qc/{sample}.failGenes"
	params:
		rulename	= "failedExons",
		threshold	= lambda wildcards: config['failed_exon_params'][config['sample_captures'][wildcards.sample]][config['sample_type'][wildcards.sample]],
		tool		= NGS_PIPELINE + "/scripts/failed_Exon_Final.pl",
		batch		= config[config['host']]['job_failedExon']
	shell:	"""
	#######################
	perl {params.tool} {input.depth_file} {params.threshold} {output.failed_exon_stats} {output.failed_gene_stats}
	#######################
	"""
