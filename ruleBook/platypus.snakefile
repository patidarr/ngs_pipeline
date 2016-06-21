rule PLATYPUS:
	input:
		bam="{subject}/{TIME}/{sample}/{sample}.bwa.final.bam",
		ref=config["reference"],
		dbsnp=config["dbsnp"],
		interval=lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]].replace("target","targetbp"),
	output:
		vcf="{subject}/{TIME}/{sample}/calls/{sample}.Platypus.raw.vcf"
	version: config["platypus"]
	log: "log/platypus.{subject}"
	params:
		rulename = "PLAT",
		batch    = config[config['host']]["job_platypus"]
	shell: """
	#######################
	module load platypus/{version}
	platypus callVariants --nCPU=${{THREADS}} --bufferSize=1000000 --maxReads=100000000 --bamFiles={input.bam} --regions={input.interval} --output={output.vcf} --refFile={input.ref}  --logFileName={log}
	sed -i 's/.bwa.final//g' {output.vcf}
	#######################
	"""
