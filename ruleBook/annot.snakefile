#ANNO_FILES =["gene","clinseq","cadd","cosmic","pcg","clinvar","hgmd","match","mcg","docm","candl","tcc","civic","sift.out","coding.rare.txt","annotations.final.txt"]
#ALL_QC = expand("{subject}/{TIME}/annotation/AnnotationInput.{suff}", suff=ANNO_FILES, TIME=TIME, subject=PATIENTS)
############
#	table_annovar for:
#		refGene
#		cytoBand
#		snp138
#		1000g2014oct_all
#		1000g2014oct_eur
#		1000g2014oct_afr
#		1000g2014oct_amr
#		1000g2014oct_eas
#		1000g2014oct_sas
#		esp6500_all
#		esp6500_ea
#		esp6500_aa
#		exac03nontcga
#		exac03
#		cg69
#		nci60
############
rule Annovar_Gene:
	input:
		config["annovar_data"]+config["annot_refgene"],
		config["annovar_data"]+config["annot_mrna"],
		config["annovar_data"]+config["annot_cb"],
		config["annovar_data"]+config["annot_snp138"],
		config["annovar_data"]+config["annot_all"],
		config["annovar_data"]+config["annot_eur"],
		config["annovar_data"]+config["annot_afr"],
		config["annovar_data"]+config["annot_amr"],
		config["annovar_data"]+config["annot_eas"],
		config["annovar_data"]+config["annot_sas"],
		config["annovar_data"]+config["annot_espall"],
		config["annovar_data"]+config["annot_espea"],
		config["annovar_data"]+config["annot_espaa"],
		config["annovar_data"]+config["annot_exacnon"],
		config["annovar_data"]+config["annot_exac"],
		config["annovar_data"]+config["annot_cg69"],
		config["annovar_data"]+config["annot_nci60"],
		file="{base}/AnnotationInput.anno",
	output:
		"{base}/AnnotationInput.gene"
	version: config["annovar"]
	params:
		rulename   = "Annot_gene",
		batch      = config[config['host']]["job_annovar"],
		RefData    = config["annovar_data"],
		build      = config["build"],
	shell: """
	#######################
	module load annovar/{version}
	table_annovar.pl {input.file} {params.RefData} -buildver {params.build} -out {input.file}.gene -remove -protocol refGene,cytoBand,snp138,1000g2014oct_all,1000g2014oct_eur,1000g2014oct_afr,1000g2014oct_amr,1000g2014oct_eas,1000g2014oct_sas,esp6500_all,esp6500_ea,esp6500_aa,exac03nontcga,exac03,cg69,nci60 -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring "-1"
	mv {input.file}.gene.{params.build}_multianno.txt {output}
	sed -i '1s/\./_/g' {output}
	rm -rf {input.file}.gene.refGene.invalid_input {input.file}.gene.invalid_input
	#######################
	"""
############
# annotate_variation
#	custom ref input. clinseq
############
rule Annovar_Clinseq:
	input:
		file="{base}/AnnotationInput.anno",
		clinseq=config["annovar_data"]+config["annot_clinseq"]
	output: "{base}/AnnotationInput.clinseq"
	version: config["annovar"]
	params:
		rulename   = "Annot_clinseq",
		batch      = config[config['host']]["job_annot"],
		RefData    = config["annovar_data"],
		build      = config["build"],
	shell: """
	#######################
	module load annovar/{version}
	annotate_variation.pl {input.file} {params.RefData} -buildver {params.build} -otherinfo --outfile {input.file}.clinseq -filter -dbtype generic -genericdbfile `basename {input.clinseq}`
	awk '{{OFS="\\t"}};{{print $3,$4,$5,$6,$7,$2}}' {input.file}.clinseq.{params.build}_generic_dropped |sed -e 's/,/\\t/g' >{output}
	head -1 {input.clinseq} >>{output}
	rm -rf {input.file}.clinseq.{params.build}_generic_dropped {input.file}.clinseq.{params.build}_generic_filtered {input.file}.clinseq.invalid_input {input.file}.clinseq.log
	#######################
	"""
################################
# Add CADD annotation
#
################################
rule Annovar_CADD:
	input:
		file="{base}/AnnotationInput.anno",
		cadd=config["annovar_data"]+config["annot_cadd"],
		cadd_indel=config["annovar_data"]+config["annot_caddind"]
	output:
		"{base}/AnnotationInput.cadd",
	version: config["annovar"]
	params:
		rulename   = "Annot_CADD",
		batch      = config[config['host']]["job_annovar"],
		RefData    = config["annovar_data"],
		build      = config["build"],
	shell: """
	#######################
	module load annovar/{version}
	annotate_variation.pl {input.file} {params.RefData} -buildver {params.build} -otherinfo -filter -dbtype cadd
	annotate_variation.pl {input.file}.{params.build}_cadd_filtered {params.RefData} -buildver {params.build} -otherinfo -filter -dbtype caddindel
	cut -f 2-7 {input.file}.{params.build}_cadd_dropped {input.file}.{params.build}_cadd_filtered.{params.build}_caddindel_dropped |sed -e 's/,/\\t/g' |awk '{{OFS="\\t"}};{{print $3,$4,$5,$6,$7,$1,$2}}' >{output}
	head -1 {input.cadd_indel} >>{output}
	rm -rf {input.file}.{params.build}_cadd_dropped {input.file}.{params.build}_cadd_filtered {input.file}.{params.build}_cadd_filtered.log {input.file}.{params.build}_cadd_filtered.{params.build}_caddindel_filtered {input.file}.{params.build}_cadd_filtered.{params.build}_caddindel_dropped {input.file}.invalid_input
	
	#######################
	"""
################################
# Add COSMIC
#
################################
rule Annovar_COSMIC:
	input:
		file="{base}/AnnotationInput.anno",
		cosmic=config["annovar_data"]+config["annot_cosmic"]
	output:
		"{base}/AnnotationInput.cosmic",
	version: config["annovar"]
	params:
		rulename   = "Annot_COSMIC",
		batch      = config[config['host']]["job_annot"],
		RefData    = config["annovar_data"],
		build      = config["build"],
	shell: """
	#######################
	module load annovar/{version}
	table_annovar.pl {input.file} {params.RefData} -buildver {params.build} --dot2underline -out {input.file}.cosmic -remove -protocol cosmic76 -operation f -nastring "NA" 
	mv {input.file}.cosmic.{params.build}_multianno.txt {output}
	rm -rf {input.file}.cosmic.invalid_input
	#######################
	"""
################################
# Add PCG
#
################################
rule Annovar_PCG:
	input:
		file="{base}/AnnotationInput.anno",
		pcg=config["annovar_data"]+config["annot_pcg"]
	output:
		"{base}/AnnotationInput.pcg",
	version: config["annovar"]
	params:
		rulename   = "Annot_PCG",
		batch      = config[config['host']]["job_annot"],
		RefData    = config["annovar_data"],
		build      = config["build"],
	shell: """
	#######################
	module load annovar/{version}
	annotate_variation.pl {input.file} {params.RefData} -buildver {params.build} -otherinfo --outfile {input.file}.pcg -filter -dbtype generic -genericdbfile `basename {input.pcg}`
	awk -F "\\t" '{{OFS="\\t"}};{{print $3,$4,$5,$6,$7,$2}}' {input.file}.pcg.{params.build}_generic_dropped |sed -e 's/,/\\t/g' >{output}
	head -1 {input.pcg} >>{output}
	rm -rf {input.file}.pcg.{params.build}_generic_dropped {input.file}.pcg.{params.build}_generic_filtered {input.file}.pcg.invalid_input {input.file}.pcg.log
	#######################
	"""
################################
# Add HGMD
#
################################
rule Annot_Custom:
	input:
		tool	=NGS_PIPELINE + "/scripts/addAnnotation.pl",
		file	="{base}/AnnotationInput.anno",
		clinvar	=config["annovar_data"]+config["annot_clinvar"],
		hgmd	=config["annovar_data"]+config["annot_hgmd"],
		match	=config["annovar_data"]+config["annot_match"],
		mcg	=config["annovar_data"]+config["annot_mcg"],
		docm	=config["annovar_data"]+config["annot_docm"],
		candl	=config["annovar_data"]+config["annot_candl"],
		tcc	=config["annovar_data"]+config["annot_tcc"],
		civic	=config["annovar_data"]+config["annot_civic"]
	output:
		clinvar="{base}/AnnotationInput.clinvar",
		hgmd   ="{base}/AnnotationInput.hgmd",
		match  ="{base}/AnnotationInput.match",
		mcg    ="{base}/AnnotationInput.mcg",
		docm   ="{base}/AnnotationInput.docm",
		candl  ="{base}/AnnotationInput.candl",
		tcc    ="{base}/AnnotationInput.tcc",
		civic  ="{base}/AnnotationInput.civic"	
	version: config["annovar"]
	params:
		rulename   = "Annot_PCG",
		batch      = config[config['host']]["job_annot"],
		RefData    = config["annovar_data"],
		build      = config["build"],
	shell: """
	#######################
	module load annovar/{version}
	{input.tool} {input.clinvar}	{input.file} >{output.clinvar}
	{input.tool} {input.hgmd}	{input.file} >{output.hgmd}
	{input.tool} {input.match}	{input.file} >{output.match}
	{input.tool} {input.mcg}	{input.file} >{output.mcg}
	{input.tool} {input.docm}	{input.file} >{output.docm}
	{input.tool} {input.candl}	{input.file} >{output.candl}
	{input.tool} {input.tcc}	{input.file} >{output.tcc}
	{input.tool} {input.civic}	{input.file} >{output.civic}	
	#######################
	"""
############
#       SIFT
############
rule SIFT:
	input:
		sift	="{base}/AnnotationInput.sift",
		tool	=NGS_PIPELINE + "/scripts/ParseSIFT.pl"
	output:
		"{base}/AnnotationInput.sift.out"
	version: config["SIFT"]
	resources: SIFT=1
	params:
		rulename   = "SIFT",
		batch      = config[config['host']]["job_SIFT"],
		build      = config["SIFTbuild"]
	shell: """
	#######################
	if [ -s {input.sift} ]; then
		module load python/2.7.9
		module load SIFT/{version}
		DIR=`pwd`
		cd ${{DIR}}/`dirname {input.sift}`
		FILE=`basename {input.sift}`
		#LOCAL=/projects/scratch/
		SIFT_exome_nssnvs.pl -i ${{FILE}} -d $SIFTDB/Human_db_37 -o ${{LOCAL}}  -z ${{DIR}}/`dirname {input.sift}`/${{FILE}}.sift_predictions.tsv
		perl {input.tool} ${{DIR}}/`dirname {input.sift}`/${{FILE}}.sift_predictions.tsv >${{DIR}}/`dirname {input.sift}`/${{FILE}}.out
	else
		echo -e "Chr\\tStart\\tEnd\\tRef\\tAlt\\tSIFT Prediction\\tSIFT Score" >{output}
	fi
	rm -rf ${{DIR}}/{input.sift}.sift_predictions.tsv
	#######################
	"""
############
#	Combine Annotation
############
rule CombineAnnotation:
	input:
		"{base}/AnnotationInput.gene",
		"{base}/AnnotationInput.clinseq",
		"{base}/AnnotationInput.cadd",
		"{base}/AnnotationInput.cosmic",
		"{base}/AnnotationInput.pcg",
		"{base}/AnnotationInput.clinvar",
		"{base}/AnnotationInput.hgmd",
		"{base}/AnnotationInput.match",
		"{base}/AnnotationInput.mcg",
		"{base}/AnnotationInput.docm",
		"{base}/AnnotationInput.candl",
		"{base}/AnnotationInput.tcc",
		"{base}/AnnotationInput.civic",
		"{base}/AnnotationInput.sift.out",
		convertor   = NGS_PIPELINE + "/scripts/CombineAnnotations.pl",
		geneanno    = NGS_PIPELINE + "/scripts/GeneAnnotation.pl",
		filter      = NGS_PIPELINE + "/scripts/filterVariants.v1.pl",
		coding      = NGS_PIPELINE + "/scripts/ProteinCoding.pl",
		blacklisted = config["annovar_data"]+ "hg19_blacklistedSites.txt"
	output: 
		filtered="{base}/AnnotationInput.coding.rare.txt",
		all	="{base}/AnnotationInput.annotations.final.txt"
	version: "1.0"
	params:
		rulename   = "combine",
		batch	   = config[config['host']]["job_Combine"],
		dataDir    = config["annovar_data"]
	shell: """
	#######################
echo "{wildcards.base}/AnnotationInput
{wildcards.base}/AnnotationInput.gene
{wildcards.base}/AnnotationInput.clinseq
{wildcards.base}/AnnotationInput.cadd
{wildcards.base}/AnnotationInput.sift.out
{wildcards.base}/AnnotationInput.clinvar
{wildcards.base}/AnnotationInput.cosmic
{wildcards.base}/AnnotationInput.hgmd
{wildcards.base}/AnnotationInput.match
{wildcards.base}/AnnotationInput.docm
{wildcards.base}/AnnotationInput.candl
{wildcards.base}/AnnotationInput.tcc
{wildcards.base}/AnnotationInput.mcg
{wildcards.base}/AnnotationInput.civic
{wildcards.base}/AnnotationInput.pcg" >{wildcards.base}/list
	perl {input.convertor} {wildcards.base}/list >{output.all}.tmp
	perl {input.geneanno} {params.dataDir}hg19_ACMG.txt {output.all}.tmp >{output.all}

	perl {input.coding} {wildcards.base}/AnnotationInput.annotations.final.txt | perl {input.filter} - {input.blacklisted} 0.05 >{output.all}.tmp
	grep -P "Chr\\tStart\\tEnd\\tRef\\tAlt" {output.all}.tmp >{output.filtered}
	grep -v -P "Chr\\tStart\\tEnd\\tRef\\tAlt" {output.all}.tmp >>{output.filtered}
	rm -rf {output.all}.tmp {wildcards.base}/list
	#######################
	"""
