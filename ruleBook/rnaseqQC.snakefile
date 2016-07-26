rule RNASeqQC:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.tophat.final.bam",
		bai="{base}/{TIME}/{sample}/{sample}.tophat.final.bam.bai",
		rna_interval=config['rRNA_interval'],
		ref_flat=config['ref_flat']
	output:
		table="{base}/{TIME}/{sample}/qc/{sample}.RnaSeqMetrics.txt",
		pdf="{base}/{TIME}/{sample}/qc/{sample}.RnaSeqMetrics.pdf"
	version: config["picard"]
	params:
		rulename  = "RnaSeqMetrics",
		batch     = config[config['host']]["job_markdup"],
	shell: """
	#######################
	module load picard/{version}
	module load R
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $PICARD_JAR CollectRnaSeqMetrics STRAND_SPECIFICITY=NONE VALIDATION_STRINGENCY=SILENT REF_FLAT={input.ref_flat} RIBOSOMAL_INTERVALS={input.rna_interval} INPUT={input.bam} OUTPUT={output.table} CHART_OUTPUT={output.pdf}
	#######################
	"""
