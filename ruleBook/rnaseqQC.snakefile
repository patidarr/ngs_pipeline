rule RNASeqQC:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.tophat.final.bam",
		bai="{base}/{TIME}/{sample}/{sample}.tophat.final.bam.bai",
		rna_interval=config['rRNA_interval'],
		ref_flat=config['ref_flat'],
	output:
		table="{base}/{TIME}/{sample}/qc/{sample}.RnaSeqMetrics.txt",
		pdf="{base}/{TIME}/{sample}/qc/{sample}.RnaSeqMetrics.pdf",
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

rule RNASeqQC1:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.tophat.final.bam",
		bai="{base}/{TIME}/{sample}/{sample}.tophat.final.bam.bai",
	output:
		table="{base}/{TIME}/{sample}/qc/{sample}.AlignmentSummaryMetrics.txt",
	version: config["picard"]
	params:
		rulename  ="RnaSeqMetrics",
		batch     =config[config['host']]["job_markdup"],
		ref       =config['Bowtie2Index'].replace('genome', 'genome.fa')
	shell: """
	#######################
	module load picard/{version}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $PICARD_JAR CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE={params.ref} INPUT={input.bam} OUTPUT={output.table} ADAPTER_SEQUENCE=null
	#######################
	"""


rule RNASeqQC_1:
	input:
		"{base}/{TIME}/{sample}/qc/fastqc/{sample}_R1_fastqc.html",
		file1="{base}/{TIME}/{sample}/qc/{sample}.RnaSeqMetrics.txt",
		file2="{base}/{TIME}/{sample}/qc/{sample}.AlignmentSummaryMetrics.txt",
		convertor=NGS_PIPELINE+ "/scripts/rnaseqQC.pl"	
	output:
		"{base}/{TIME}/{sample}/qc/{sample}.RnaSeqQC.txt"
	params:
		rulename  = "RnaSeqMetrics",
		batch     = config[config['host']]["job_default"],
		diagnosis = lambda wildcards: config['Diagnosis'][wildcards.sample]
	shell: """
	#######################
	perl {input.convertor} {wildcards.base}/{TIME}/{wildcards.sample}/qc/fastqc/{wildcards.sample}_R1_fastqc/fastqc_data.txt {input.file2} {input.file1} {wildcards.base} {wildcards.sample} "{params.diagnosis}" >{output}
	#######################
	"""
rule RNASeqQC_2:
	input : lambda wildcards: SUB_QC[wildcards.subject]
	output: "{subject}/{TIME}/qc/{subject}.RnaSeqQC.txt"
	params:
		rulename  = "QC_Sum",
		batch     = config[config['host']]["job_default"]
	shell: """
	#######################
	export LC_ALL=C
        cat {input} |sort |uniq |sed -e '/^$/d'>{output}	
	#######################
	"""
rule RNASeqQC_3:
	input : RNA_QC_ALL
	output: "RnaSeqQC.txt"
	params:
		rulename  = "QC_Sum",
		batch     = config[config['host']]["job_default"]
	shell: """
	#######################
        export LC_ALL=C
        cat {input} {output} |sort|uniq |sed -e '/^$/d'>{output}
	#######################
	"""
