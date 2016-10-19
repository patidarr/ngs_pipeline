############
#       GenotypeFile
############
# Using Older version of samtools for this purpose
rule PILE4SEQ:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.bwa.final.bam",
		ref=config["reference"],
		seq=NGS_PIPELINE+ "/scripts/sequenza-utils.py"
	output: temp("{base}/{TIME}/{sample}/{sample}.mpileup.gz")
	version: config["samtools_old"]
	params:
		rulename  = "pile4seq",
		batch     = config[config['host']]["job_genotype"],
	shell: """
	#######################
	module load samtools/{version}
	samtools mpileup -Q 20 -q 30 -f {input.ref} {input.bam}| gzip > {output}
	#######################
	"""
rule Sequenza:
	input:
		files=lambda wildcards: SequenzaPairs[wildcards.Tumor],
		seq=NGS_PIPELINE+ "/scripts/sequenza-utils.py",
		gc_ref=config["annovar_data"]+config["gc50Base"],
		RCode=NGS_PIPELINE+ "/scripts/run_sequenza_pipeline.R",
	output:
		all=temp("{subject}/{TIME}/{Tumor}/sequenza/{Tumor}.seqz.gz"),
		bin="{subject}/{TIME}/{Tumor}/sequenza/{Tumor}.seqz_small.gz"
	version: config['pypy']
	params:
		rulename = "sequenza",
		batch    = config[config['host']]['job_genotype']
	shell: """
	#######################
	module load pypy/{version}
	pypy {input.seq} pileup2seqz -gc {input.gc_ref} gc_ref -n {input.files[0]} -t {input.files[1]} |gzip >{output.all}
	pypy {input.seq} seqz-binning -w 50 -s {output.all} | gzip > {output.bin}
	module load R
	cd {wildcards.subject}/{wildcards.TIME}/{wildcards.Tumor}/sequenza/
	{input.RCode} {wildcards.Tumor}
	#######################
	"""
