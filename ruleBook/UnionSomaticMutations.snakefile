UNION_SOM_MUT={}
for sample in config['sample_references'].keys():
	subject=SAMPLE_TO_SUBJECT[sample]
	local =[(subject+"/"+TIME+"/"+sample+"/calls/"+sample+".strelka.snvs.annotated.txt"),
	        (subject+"/"+TIME+"/"+sample+"/calls/"+sample+".strelka.indels.annotated.txt"),
	        (subject+"/"+TIME+"/"+sample+"/calls/"+sample+".MuTect.annotated.txt")]
	if sample in config['sample_RNASeq'].keys():
		local = [w.replace('annotated','annotated.expressed') for w in local]
	UNION_SOM_MUT[sample] = local
	TARGET +=[subject+"/"+TIME+ACT_DIR+sample+".unionSomaticVars.txt"]
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
