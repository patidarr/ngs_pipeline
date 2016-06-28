############
#       Union Somatic Calls
############
rule UnionSomaticCalls:
	input:
		Files=lambda wildcards: UNION_SOM_MUT[wildcards.sample]
	output: 
		"{base}/{sample}.unionSomaticVars.txt"
	params:
		rulename = "UnionSomaticCalls",
		batch    = config[config['host']]['job_default'],
		tool     = NGS_PIPELINE + "/scripts/UnionSomaticCalls.pl"
	shell: """
	#######################
	perl {params.tool} {input.Files} >{output}
	#######################
	"""
