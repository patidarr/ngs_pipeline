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
		batch    = config[config['host']]['job_Sequenza']
	shell: """
	#######################
	module load pypy/{version}
	pypy {input.seq} pileup2seqz -gc {input.gc_ref} -n {input.files[0]} -t {input.files[1]} |gzip >{output.all}
	pypy {input.seq} seqz-binning -w 50 -s {output.all} | gzip > {output.bin}
	module load R
	cd {wildcards.subject}/{wildcards.TIME}/{wildcards.Tumor}/sequenza/
	{input.RCode} --sample {wildcards.Tumor}
	#######################
	"""

rule Sequenza_geneAnnot:
	input:
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}/{Tumor}_segments.txt"
	output:
		"{subject}/{TIME}/{Tumor}/sequenza/{Tumor}.txt"
	params:
		rulename = "sequenza_1",
		batch    = config[config['host']]['job_Sequenza']
	shell: """
	#######################
	#module load bedtools/2.25.0
	##grep -v "chromosome" OM161/20160415/CL0034_T_E/sequenza/CL0034_T_E/CL0034_T_E_segments.txt |sed -e 's/"//g'|intersectBed -a /projects/Clinomics/Tools/ref/hg19/gene_coordinates.txt -b -
	##sed -e 's/"//g' OM161/20160415/CL0034_T_E/sequenza/CL0034_T_E/CL0034_T_E_segments.txt |/projects/Clinomics/Tools/ngs_pipeline//scripts/AddGene.pl out.tmp -
		
	#######################
	"""
