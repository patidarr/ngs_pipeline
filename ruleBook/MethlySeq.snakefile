if 'methylseq' in config:
	for subject in config['methylseq'].keys():
		if subject not in PATIENTS:
			PATIENTS.append(subject)
		for library in config['methylseq'][subject]:
			ALL_FASTQC += [subject+"/"+TIME+"/"+library+"/qc/fastqc/"+library+"_R2_fastqc.html"]
			ALL_QC += [subject+"/"+TIME+"/"+library+"/"+library+"_pe.deduplicated.CpG_report.txt.gz"]
			ALL_QC += [subject+"/"+TIME+"/"+library+"/"+library+"_PE_report.html"]
	for subject in config['methylseq']:
		SUBJECT_TO_SAMPLE[subject] = expand("{sample}", sample = config['methylseq'][subject])
############
#       Bismark
############
rule Bismark:
	input: 
		R=lambda wildcards: FQ[wildcards.sample],
		ref=config["BismarkIndex"]
	output:
		"{base}/{TIME}/{sample}/{sample}_pe.bam",
	version: config["bismark"]
	params:
		rulename  = "bismark",
		platform  = config["platform"],
		samtools  = config["samtools"],
		batch     = config[config['host']]["job_bismark"]
	shell: """
	#######################
	module load bismark/{version} bowtie
	module load samtools/{params.samtools}
	bismark -p ${{THREADS}} {input.ref} -1 {input.R[0]} -2 {input.R[1]} --output_dir {wildcards.base}/{TIME}/{wildcards.sample}/ --basename {wildcards.sample} --temp_dir ${{LOCAL}}
	#######################
	"""
############
#       Bismark Deduplicate
############
rule BisDedup:
	input:
		"{base}/{TIME}/{sample}/{sample}_pe.bam"
	output:
		"{base}/{TIME}/{sample}/{sample}_pe.deduplicated.bam"
	version: config["bismark"]
	params:
		rulename  = "bisdedup",
		samtools  = config["samtools"],
		batch     = config[config['host']]["job_bismark"]
	shell: """
	#######################
	module load bismark/{version} bowtie samtools/{params.samtools}
	deduplicate_bismark -p {input} --bam
	#######################
	"""
############
#       Bismark Methylation Extractor
############
rule BismarkMethExt:
	input:
		bam="{base}/{TIME}/{sample}/{sample}_pe.deduplicated.bam",
		ref=config["BismarkIndex"]
	output:
		"{base}/{TIME}/{sample}/{sample}_pe.deduplicated.CpG_report.txt.gz",
	version: config["bismark"]
	params:
		rulename  = "BisMetExt",
		samtools  = config["samtools"],
		batch     = config[config['host']]["job_bismark"]
	shell: """
	#######################
	module load bismark/{version} bowtie
	module load samtools/{params.samtools}
	bismark_methylation_extractor --gzip {input.bam} --ample_memory  --multicore ${{THREADS}} --output {wildcards.base}/{TIME}/{wildcards.sample}/ --genome_folder {input.ref} --bedGraph --cytosine_report --report 
	#bismark_methylation_extractor {input.bam} --genome_folder {input.ref} --ample_memory  --multicore ${{THREADS}} --report --merge_non_CpG --output {wildcards.base}/{TIME}/{wildcards.sample}/ --comprehensive 
	#######################
	"""
############
#	Bismark Nucleotide Coverage
############
rule bam2nuc:
	input: "{base}/{TIME}/{sample}/{sample}_pe.bam"
	output: "{base}/{TIME}/{sample}/{sample}.nucleotide_stats.txt"
	version: config["bismark"]
	params:
		rulename  = "BisBam2Nuc",
		samtools  = config["samtools"],
		batch     = config[config['host']]["job_bismark"],
		ref=config["BismarkIndex"]
	shell: """
	#######################
	module load bismark/{version} samtools/{params.samtools}
	bam2nuc --dir {wildcards.base}/{TIME}/{wildcards.sample}/ --genome_folder {params.ref} {input}
	#######################
	"""
############
#       Bismark HTML Report
############
rule BisReport: 
	input:
		"{base}/{TIME}/{sample}/{sample}_pe.deduplicated.CpG_report.txt.gz",
		"{base}/{TIME}/{sample}/{sample}.nucleotide_stats.txt"
	output:
		"{base}/{TIME}/{sample}/{sample}_PE_report.html"
	version: config["bismark"]
	params:
		rulename  = "BisReport",
		samtools  = config["samtools"],
		batch     = config[config['host']]["job_default"]
	shell: """
	#######################
	module load bismark/{version} samtools/{params.samtools}
	bismark2report --dir {wildcards.base}/{TIME}/{wildcards.sample}/ --alignment_report {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}_PE_report.txt --dedup_report {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}_pe.deduplication_report.txt --splitting_report {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}_pe.deduplicated_splitting_report.txt --mbias_report {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}_pe.deduplicated.M-bias.txt
	#bismark2report --dir {wildcards.base}/{TIME}/{wildcards.sample}/ --alignment_report {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}_PE_report.txt --dedup_report {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}_pe.deduplication_report.txt --splitting_report {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}_pe.deduplicated_splitting_report.txt --mbias_report {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}_pe.deduplicated.M-bias.txt --nucleotide_report {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}
	#######################
	"""
