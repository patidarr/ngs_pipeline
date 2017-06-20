TARGET     += ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.squeeze.done".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
for subject  in config['RNASeq'].keys():
	for sample in config['RNASeq'][subject]:
		TARGET      += [subject+"/"+TIME+"/"+sample+"/qc/"+sample+".star.squeeze.done"]
############
#       bamUtil
############
rule bamUtil:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.{aligner}.final.bam",
		ref=config["reference"],
	output: "{base}/{TIME}/{sample}/qc/{sample}.{aligner}.squeeze.done"
	version: config["bamutil"]
	params:
		rulename  = "bamutil",
		samtools  = config['samtools'],
		batch     = config[config['host']]["job_bamUtil"]
	shell: """
	#######################
	module load bamutil/{version}
	module load samtools/{params.samtools}
	bam squeeze --in {input.bam} --out {wildcards.base}/{wildcards.TIME}/{wildcards.sample}/{wildcards.sample}.{wildcards.aligner}.final.squeeze.bam --refFile {input.ref} --rmTags "PG:Z;RG:Z;BI:Z;BD:Z"
	samtools index {wildcards.base}/{wildcards.TIME}/{wildcards.sample}/{wildcards.sample}.{wildcards.aligner}.final.squeeze.bam
	touch {output}
	#######################
        """
