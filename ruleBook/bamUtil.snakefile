
############
#       bamUtil
############
rule bamUtil:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.{aligner}.final.bam",
		ref=config["reference"],
	output: "{base}/{TIME}/{sample}/{sample}.{aligner}.final.squeeze.bam"
	version: config["bamutil"]
	params:
		rulename  = "bamutil",
		batch     = config[config['host']]["job_annot"]
	shell: """
	#######################
	module load bamutil/{version}
	module load samtools
	bam squeeze --in {input.bam} --out {output} --refFile {input.ref} --rmTags "PG:Z;RG:Z;BI:Z;BD:Z"
	samtools index {output}
	#######################
        """
