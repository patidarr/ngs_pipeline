rule TargetIntervals:
	input:
		bam="{base}/{sample}/{sample}.bwa.final.bam",
		bai="{base}/{sample}/{sample}.bwa.final.bam.bai"
	output:
		probe_intervals  = temp("{base}/qc/{sample}.probe.intervals"),
		target_intervals = temp("{base}/qc/{sample}.target.intervals")
	version:
		config['samtools']
	params:
		rulename = "targetIntervals",
		target_intervals=lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]],
		probe_intervals=lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]].replace(".target.", ".design."),
		batch	= config[config['host']]["job_default"],
	shell:	"""
	#######################
	module load samtools/{version}
	cat <(samtools view -H {input.bam}) <(gawk '{{print $1 "\t" $2+1 "\t" $3 "\t+\tinterval_" NR}}' {params.probe_intervals}  )> {output.probe_intervals}
	cat <(samtools view -H {input.bam}) <(gawk '{{print $1 "\t" $2+1 "\t" $3 "\t+\tinterval_" NR}}' {params.target_intervals} )> {output.target_intervals} 
	#######################
	"""

rule HSMetrics:
	input:
		bam="{base}/{sample}/{sample}.bwa.final.bam",
		probe_intervals  = "{base}/qc/{sample}.probe.intervals",
		target_intervals = "{base}/qc/{sample}.target.intervals",
	output:
		"{base}/qc/{sample}.hsmetrics"
	version: config['picard']
	params:
		rulename = "hsMetrics",
		reference = config["reference"],
		batch	= config[config['host']]["job_markdup"],
	shell: """
	#######################
	module load picard/{version}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $PICARD_JARPATH/CalculateHsMetrics.jar BAIT_INTERVALS={input.probe_intervals} TARGET_INTERVALS={input.target_intervals} INPUT={input.bam} OUTPUT={output} METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE={params.reference} QUIET=true  VALIDATION_STRINGENCY=SILENT
	#######################
	"""
