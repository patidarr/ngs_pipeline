############
#	Mutational Signature
############
rule MutationalSignature:
	input:
		File ="{base}/{sample}.unionSomaticVars.txt",
		tool   =NGS_PIPELINE + "/scripts/mutationSignature.R"
	output: 
		v1="{base}/{sample}.mutationalSignature.pdf"
	version: config["R"]
	params:
		rulename = "MutationalSignature",
		batch    = config[config['host']]["job_default"],
	shell: """
	#######################
	module load R
	awk '{{OFS="\\t"}}{{print $1,$2,$4,$5,"{wildcards.sample}"}}' {input.File} >{output.v1}.tmp
	{input.tool} --input {output.v1}.tmp --sample {wildcards.sample} --output {output.v1}
	rm -rf {output.v1}.tmp
	#######################
	"""
