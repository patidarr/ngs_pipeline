TARGET     += ["{subject}/{TIME}/{sample}/verifyBamID/{sample}.selfSM".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
########################
#	verifyBamID 1.1.3 
########################
rule verifyBamID:
	input:
		bam="{subject}/{TIME}/{sample}/{sample}.bwa.final.bam",
		vcf=lambda wildcards: config['verifyBamID_vcf'][config['sample_captures'][wildcards.sample]],
	output:
		"{subject}/{TIME}/{sample}/verifyBamID/{sample}.selfSM"
	version: config["verifybamid"]
	params:
		rulename  = "verifyBamID",
		batch	  = config[config['host']]["job_default"]
	shell: """
	#######################
	module load verifybamid/{version}
	verifyBamID --vcf {input.vcf} --bam {input.bam} --maxDepth 3000 --ignoreRG --site --chip-none --precise --minMapQ 30 --minQ 20 --minAF 0.05 --out {wildcards.subject}/{wildcards.TIME}/{wildcards.sample}/verifyBamID/{wildcards.sample}
	#######################
	"""
