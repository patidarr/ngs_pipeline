############
#	pVACSeq
############
rule pVACSeq:
	input:
		files=lambda wildcards: UNION_SOM_MUT[wildcards.sample],
		tool =NGS_PIPELINE + "/scripts/consensusSomaticVCF.pl"
	output: 
		vcf="{base}/{TIME}/{sample}/HLA/{sample}.somatic.vep.vcf"
	version: config["R"]
	params:
		rulename = "pVACSeq",
		normal	 = lambda wildcards: config['sample_references'][wildcards.sample][0],
		batch    = config[config['host']]["job_default"]
	shell: """
	#######################
	module load vcftools vep
	perl {input.tool} -vcf {wildcards.base}/{wildcards.TIME}/{wildcards.sample}/calls/{wildcards.sample}.strelka.indels.raw.vcf,{wildcards.base}/{wildcards.TIME}/{wildcards.sample}/calls/{wildcards.sample}.strelka.snvs.raw.vcf,{wildcards.base}/{wildcards.TIME}/{wildcards.sample}/calls/{wildcards.sample}.MuTect.raw.vcf -order {params.normal},{wildcards.sample} -filter REJECT |vcf-subset -u -c {wildcards.sample} >{output.vcf}
	perl /apps/VEP/81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i CL0062_T1D_EM.MuTect.raw.vcf --plugin Downstream --plugin Wildtype --terms SO --offline --cache --dir_cache /apps/VEP/81/cache/ --assembly GRCh37 --output_file CL0062_T1D_EM.MuTect.raw.vep.vcf --vcf --force_overwrite
	#######################
	"""
