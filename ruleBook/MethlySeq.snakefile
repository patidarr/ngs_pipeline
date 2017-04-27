if 'methylseq' in config:
	for subject in config['methylseq'].keys():
		if subject not in PATIENTS:
			PATIENTS.append(subject)
		for library in config['methylseq'][subject]:
			ALL_FASTQC += [subject+"/"+TIME+"/"+library+"/qc/fastqc/"+library+"_R2_fastqc.html"]
			ALL_QC += [subject+"/"+TIME+"/"+library+"/"+library+".bismark.bam"]
			ALL_QC += [subject+"/"+TIME+"/"+library+"/"+library+".bismark.report.txt"]
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
		"{base}/{TIME}/{sample}/{sample}.bismark.bam",
		"{base}/{TIME}/{sample}/{sample}.bismark.bam.bai"
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
	mv {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}_pe.bam {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}.bismark.bam
	samtools index {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}.bismark.bam
	#######################
	"""
############
#       Bismark Methylation Extractor
############
rule BismarkMethExt:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.bismark.bam",
		ref=config["BismarkIndex"]
	output:
		"{base}/{TIME}/{sample}/{sample}.bismark.report.txt",
	version: config["bismark"]
	params:
		rulename  = "BisMetExt",
		samtools  = config["samtools"],
		batch     = config[config['host']]["job_bismark"]
	shell: """
	#######################
	module load bismark/{version} bowtie
	module load samtools/{params.samtools}
	bismark_methylation_extractor {input.bam} --genome_folder {input.ref} --ample_memory  --multicore ${{THREADS}} --report --merge_non_CpG --output {wildcards.base}/{TIME}/{wildcards.sample}/ --comprehensive 
	#######################
	"""
