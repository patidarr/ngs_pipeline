SequenzaPairs ={}
if len(config['sample_references']) > 0:
	for Tumor in config['sample_references']:
		for Normal in config['sample_references'][Tumor]:
			if config['sample_captures'][Tumor] not in config['Panel_List']:
				SequenzaPairs[Tumor] = ["{subject}/{TIME}/{sample}/{sample}.mpileup.gz".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[Normal], sample=Normal), "{subject}/{TIME}/{sample}/{sample}.mpileup.gz".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[Tumor], sample=Tumor) ]
for sample in config['sample_references'].keys():
	subject=SAMPLE_TO_SUBJECT[sample]
	if config['sample_captures'][sample] not in config['Panel_List']:
		TARGET +=[subject+"/"+TIME+"/"+sample+"/sequenza/"+sample+"/"+sample+"_alternative_fit.pdf"]	
		TARGET +=[subject+"/"+TIME+"/"+sample+"/sequenza/"+sample+".txt"]

############
#       GenotypeFile
############
# Using Older version of samtools for this purpose
rule PILE4SEQ:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.bwa.final.bam",
		ref=config["reference"],
		genome=config["reference"].replace(".fasta",".genome"),
		interval= lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]],
		seq=NGS_PIPELINE+ "/scripts/sequenza-utils.py"
	output: temp("{base}/{TIME}/{sample}/{sample}.mpileup.gz")
	version: config["samtools_old"]
	params:
		rulename  = "pile4seq",
		batch     = config[config['host']]["job_Sequenza"],
	shell: """
	#######################
	module load bedtools/2.25.0
	slopBed -i {input.interval} -g {input.genome} -b 200 |mergeBed -i - >${{LOCAL}}/Region.bed	
	module load samtools/{version}
	samtools mpileup -Q 20 -q 30 -L ${{LOCAL}}/Region.bed  -f {input.ref} {input.bam}| gzip > {output}
	#######################
	"""
rule Sequenza:
	input:
		files=lambda wildcards: SequenzaPairs[wildcards.Tumor],
		seq=NGS_PIPELINE+ "/scripts/sequenza-utils.py",
		gc_ref=config["annovar_data"]+config["gc50Base"],
		RCode=NGS_PIPELINE+ "/scripts/run_sequenza_pipeline.R",
	output:
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_CN_bars.pdf",
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_CP_contours.pdf",
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_alternative_fit.pdf",
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_alternative_solutions.txt",
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_chromosome_view.pdf",
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_confints_CP.txt",
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_genome_view.pdf",
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_model_fit.pdf",
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_mutations.txt",
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_segments.txt",
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_sequenza_cp_table.RData",
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_sequenza_extract.RData",
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_sequenza_log.txt",
		all=temp("{subject}/{TIME}/{Tumor}/sequenza/{Tumor}.seqz.gz"),
		bin=temp("{subject}/{TIME}/{Tumor}/sequenza/{Tumor}.seqz_small.gz")
	version: config['pypy']
	params:
		rulename = "sequenza",
		R	 = config["version_R"],
		batch    = config[config['host']]['job_Sequenza']
	shell: """
	#######################
	module load pypy/{version}
	pypy {input.seq} pileup2seqz -gc {input.gc_ref} -n {input.files[0]} -t {input.files[1]} |gzip >{output.all}
	pypy {input.seq} seqz-binning -w 50 -s {output.all} | gzip > {output.bin}
	module load R/{params.R}
	cd {wildcards.subject}/{wildcards.TIME}/{wildcards.Tumor}/sequenza/
	{input.RCode} --sample {wildcards.Tumor}
	#######################
	"""

rule Sequenza_geneAnnot:
	input:
		file="{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_segments.txt",
		interval=lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.Tumor]],
		convertor=NGS_PIPELINE+"/scripts/GeneAnnotation.v1.pl",
		geneList=config["annovar_data"]+config["geneList"]
	output:
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}.txt"
	version: config["bedtools"]
	params:
		rulename = "sequenza_1",
		batch    = config[config['host']]['job_Sequenza']
	shell: """
	#######################
	set +eo pipefail
        module load bedtools/{version}
        sed -e 's/"//g' {input.file} |sed -e 's/chromosome/#chromosome/' | bedtools intersect  -wa -a {input.interval} -b - |grep -v NOTFOUND |sed -e 's/___/\\t/g'| cut -f 1-4| bedtools expand -c 4 >{output}.temp
        sed -e 's/"//g' {input.file} |sed -e 's/chromosome/#chromosome/' |head -1 >{output}.temp1
        sed -i 's/end.pos\\tBf/end.pos\\tGene\\tBf/g' {output}.temp1
        sed -e 's/"//g' {input.file} |sed -e 's/chromosome/#chromosome/' |intersectBed -a {output}.temp -b - -wb |cut -f 1-4,8-100 >>{output}.temp1
        perl {input.convertor} {input.geneList} {output}.temp1 3 >{output}
        rm -rf {output}.temp {output}.temp1	
	#######################
	"""
