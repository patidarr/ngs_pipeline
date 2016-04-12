rule PLATYPUS:
	input:
		bam="{subject}/{sample}/{sample}.bwa.final.bam",
		ref=config["reference"],
		dbsnp=config["dbsnp"],
		interval=config["coding_bed"]
	output:
		vcf="{subject}/{sample}/calls/{sample}.platypus.raw.vcf"
	version: config["platypus"]
	log: "log/platypus.{subject}"
	params:
		rulename = "PLAT",
		batch    = config[config['host']]["job_platypus"]
	shell: """
	#######################
	module load platypus/{version}
	gawk '{{print $1 "\t" $2-1 "\t" $3}}' {input.interval} > ${{LOCAL}}/target_intervals.bed	
	platypus callVariants --nCPU=${{THREADS}} --bufferSize=1000000 --maxReads=100000000 --bamFiles={input.bam} --regions=${{LOCAL}}/target_intervals.bed --output={output.vcf} --refFile={input.ref}  --logFileName={log}
	sed -i 's/.bwa.final//g' {output.vcf}
	#######################
	"""
