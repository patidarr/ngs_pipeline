SUB_CON_QC = {}
for subject in config['subject'].keys():
	SUB_CON_QC[subject]  = ["{subject}/{TIME}/{sample}/qc/{sample}.consolidated_QC".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	TARGET += ["{subject}/{TIME}/qc/{subject}.consolidated_QC.txt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
CON_QC = ["{subject}/{TIME}/{sample}/qc/{sample}.consolidated_QC".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
TARGET += ["{subject}/{TIME}/{sample}/qc/{sample}.consolidated_QC".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]

############
#       QC
############
rule QC:
        input:
                bam="{base}/{TIME}/{sample}/{sample}.bwa.final.bam",
                hsMatrix="{base}/{TIME}/{sample}/qc/{sample}.hsmetrics",
		target_intervals =lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]],
                tool= NGS_PIPELINE+ "/scripts/QC_stats_Final.py",
		tool_perl       = NGS_PIPELINE+ "/scripts/addAttributes.pl",
        output: "{base}/{TIME}/{sample}/qc/{sample}.consolidated_QC"
        version:
                config['samtools']
        params:
                rulename        = "consolidated_QC",
                bedtools        = config['bedtools'],
		python		= config["version_python"],
		tool_perl       = NGS_PIPELINE+ "/scripts/addAttributes.pl",
                batch           = config[config['host']]["job_c_QC"],
                diagnosis       = lambda wildcards: config['Diagnosis'][wildcards.sample]
        shell:  """
        #######################
        module load python/{params.python}
        module load samtools/{version}
        module load bedtools/{params.bedtools}
        python {input.tool} {input.bam} {input.target_intervals} ${{LOCAL}} {wildcards.base} {wildcards.sample}  "{params.diagnosis}" > {output}.tmp
	perl   {input.tool_perl} {wildcards.sample} {input.hsMatrix} {output}.tmp  {output}
	rm -rf {output}.tmp
        #######################
        """
############
#       QC Summary for Patient
############
rule QC_Summary_Patient:
	input : lambda wildcards: SUB_CON_QC[wildcards.subject]
	output: "{subject}/{TIME}/qc/{subject}.consolidated_QC.txt"
	params:
		rulename  = "QC_Sum",
		batch     = config[config['host']]["job_default"]
	shell: """
	#######################
	export LC_ALL=C
	cat {input} |sort |uniq |sed -e '/^$/d'>{output}
	#######################
	"""
############
#	QC Summary for Cohort
############
rule QC_Summary:
	input : CON_QC
	output: "Consolidated_QC.txt"
	params:
		rulename  = "QC_Sum",
		batch	  = config[config['host']]["job_default"]
	shell: """
	#######################
	touch {output}
	export LC_ALL=C
	cat {input} {output} |sort|uniq |sed -e '/^$/d'>{output}.tmp
	mv {output}.tmp {output}
	#######################
	"""
