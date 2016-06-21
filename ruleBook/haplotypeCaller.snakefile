rule HAPCALLER:
	input: 	
		bam="{subject}/{TIME}/{sample}/{sample}.bwa.final.bam",
		bai="{subject}/{TIME}/{sample}/{sample}.bwa.final.bam.bai",
		ref=config["reference"],
		dbsnp=config["dbsnp"],
		interval=lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]].replace("target", "targetbp"),
	output:
		vcf="{subject}/{TIME}/{sample}/calls/{sample}.HC_DNASeq.raw.vcf"
	version: config["GATK"]
	params:
		rulename = "HC",
		batch    = config[config['host']]["job_gatk"]
	shell: """
	#######################
	module load GATK/{version}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.interval} -o {output.vcf} --dbsnp {input.dbsnp} -mbq 20 -mmq 30
	#######################
	"""
