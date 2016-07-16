rule QC:
        input:
                bam="{base}/{TIME}/{sample}/{sample}.bwa.final.bam",
                hsMatrix="{base}/{TIME}/{sample}/qc/{sample}.hsmetrics",
        output: "{base}/{TIME}/{sample}/qc/{sample}.consolidated_QC"
        version:
                config['samtools']
        params:
                rulename        = "consolidated_QC",
                target_intervals= lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]],
                bedtools        = config['bedtools'],
                tool            = NGS_PIPELINE+ "/scripts/QC_stats_Final.py",
		tool_perl       = NGS_PIPELINE+ "/scripts/addAttributes.pl",
                batch           = config[config['host']]["job_c_QC"],
                diagnosis       = lambda wildcards: config['Diagnosis'][wildcards.sample]
        shell:  """
        #######################
        module load python/2.7.9
        module load samtools/{version}
        module load bedtools/{params.bedtools}
        python {params.tool} {input.bam} {params.target_intervals} ${{LOCAL}} {wildcards.base} {wildcards.sample}  "{params.diagnosis}" > {output}.tmp
	perl {params.tool_perl} {wildcards.sample} {input.hsMatrix} {output}.tmp  {output}
	rm -rf {output}.tmp
        #######################
        """
rule QC_Summary_Patient:
	input : lambda wildcards: SUB_CON_QC[wildcards.subject]
	output: "{subject}/{TIME}/qc/{subject}.consolidated_QC.txt"
	params:
		rulename  = "QC_Sum",
		tools     = NGS_PIPELINE+ "/scripts/awk_sort_withHeader.awk",
		batch     = config[config['host']]["job_default"]
	shell: """
	#######################
	cat {input} |{params.tools} |uniq |sed -e '/^$/d'>{output}
	#######################
	"""
rule QC_Summary:
	input : CON_QC
	output: "Consolidated_QC.txt"
	params:
		rulename  = "QC_Sum",
		tools	  = NGS_PIPELINE+ "/scripts/awk_sort_withHeader.awk",
		batch	  = config[config['host']]["job_default"]
	shell: """
	#######################
	touch {output}
	export LC_ALL=C
	cat {input} {output} |sort|uniq |sed -e '/^$/d'>{output}.tmp
	mv {output}.tmp {output}
	#######################
	"""
