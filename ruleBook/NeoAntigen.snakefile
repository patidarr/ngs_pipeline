############
#	pVACSeq
############
rule pVACSeq:
	input:
		files=lambda wildcards: UNION_SOM_MUT[wildcards.sample],
		tool =NGS_PIPELINE + "/scripts/consensusSomaticVCF.pl"
	output: 
		vcf="{base}/{TIME}/{sample}/NeoAntigen/{sample}.somatic.vep.vcf"
	version: config["R"]
	params:
		rulename = "pVACSeq",
		normal	 = lambda wildcards: config['sample_references'][wildcards.sample][0],
		batch    = config[config['host']]["job_default"]
	shell: """
	#######################
	module load vcftools vep pvacseq
	perl {input.tool} -vcf {wildcards.base}/{wildcards.TIME}/{wildcards.sample}/calls/{wildcards.sample}.strelka.indels.raw.vcf,{wildcards.base}/{wildcards.TIME}/{wildcards.sample}/calls/{wildcards.sample}.strelka.snvs.raw.vcf,{wildcards.base}/{wildcards.TIME}/{wildcards.sample}/calls/{wildcards.sample}.MuTect.raw.vcf -order {params.normal},{wildcards.sample} -filter REJECT |vcf-subset -u -c {wildcards.sample} >{output.vcf}.tmp
	variant_effect_predictor.pl -i {output.vcf}.tmp --plugin Downstream --plugin Wildtype --terms SO --offline --cache --dir_cache $VEPCACHEDIR --assembly GRCh37 --output_file {output.vcf} --vcf --force_overwrite
	rm -rf {output.vcf}.tmp
	allele=`cut -f1 {wildcards.base}/{wildcards.TIME}/{params.normal}/HLA/{params.normal}.Calls.txt|grep -v Allele|tr '\\n' ',' |sed -e 's/,$//g'`
	pvacseq run -e 8,9,10,11 {output.vcf} {wildcards.sample} ${{allele}} {{NNalign,NetMHC,NetMHCIIpan,NetMHCcons,NetMHCpan,PickPocket,SMM,SMMPMBEC,SMMalign}} {wildcards.base}/{wildcards.TIME}/{wildcards.sample}/NeoAntigen/
	#######################
	"""
