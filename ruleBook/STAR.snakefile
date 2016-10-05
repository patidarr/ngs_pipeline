############
#       STAR
############
rule STAR_TPM:
	input:  R=lambda wildcards: FQ[wildcards.sample],
		ref=config["reference"],
		gtf1=config['GTF']['UCSC'],
		gtf2=config['GTF']['ENS']
	output:
		temp("{subject}/{TIME}/{sample}/{sample}.star_UCSC.bam"),
		temp("{subject}/{TIME}/{sample}/{sample}.star_ENS.bam")
	version: config["STAR"]
	params:
		rulename  = "STAR",
		batch     = config[config['host']]['job_STAR_TPM'],
		star_ref  = config['STAR_ref'],
		awk       = NGS_PIPELINE + "/scripts/SJDB.awk",
		home      = WORK_DIR,
	shell: """
	#######################
	module load STAR/{version}
	cd ${{LOCAL}}/
	# run 1st pass
	STAR --outTmpDir STEP1 \
		--genomeDir {params.star_ref} \
		--readFilesIn {input.R[0]} {input.R[1]} \
		--readFilesCommand zcat\
		--outSAMtype BAM SortedByCoordinate\
		--outFileNamePrefix {wildcards.sample} \
		--runThreadN ${{THREADS}} \
		--outFilterMismatchNmax 2
	echo "Finished Step 1"

	# make splice junctions database file out of SJ.out.tab, filter out non-canonical junctions
	mkdir GenomeForPass2
	awk -f {params.awk} {wildcards.sample}SJ.out.tab > GenomeForPass2/{wildcards.sample}.out.tab.Pass1.sjdb
	echo "Finished Step 2"

	# generate genome with junctions from the 1st pass
	STAR --outTmpDir STEP2\
		--genomeDir GenomeForPass2\
		--runMode genomeGenerate\
		--genomeSAindexNbases 8\
		--genomeFastaFiles {input.ref}\
		--sjdbFileChrStartEnd GenomeForPass2/{wildcards.sample}.out.tab.Pass1.sjdb\
		--sjdbOverhang 100\
		--runThreadN ${{THREADS}}
	echo "Finished Step 3"

	# run 2nd pass with the new genome
	STAR --outTmpDir STEP3\
		--genomeDir GenomeForPass2\
		--runThreadN ${{THREADS}}\
		--outSAMattributes All\
		--readFilesIn {input.R[0]} {input.R[1]}\
		--outSAMtype BAM SortedByCoordinate\
		--sjdbGTFfile {input.gtf1}\
		--readFilesCommand zcat\
		--outFileNamePrefix {wildcards.sample}_ucsc
	cp {wildcards.sample}_ucscAligned.sortedByCoord.out.bam {params.home}/{wildcards.subject}/{TIME}/{wildcards.sample}/{wildcards.sample}.star_ucsc.bam

	echo "Finished Step 4"


	STAR --outTmpDir STEP4\
		--genomeDir GenomeForPass2\
		--runThreadN ${{THREADS}}\
		--outSAMattributes All\
		--readFilesIn {input.R[0]} {input.R[1]}\
		--outSAMtype BAM SortedByCoordinate\
		--sjdbGTFfile {input.gtf2}\
		--readFilesCommand zcat\
		--outFileNamePrefix {wildcards.sample}_ens
	cp {wildcards.sample}_ensAligned.sortedByCoord.out.bam {params.home}/{wildcards.subject}/{TIME}/{wildcards.sample}/{wildcards.sample}.star_ens.bam
	echo "Finished Step 5"
	#######################
	"""
############
# featureCounts
#############
rule FeatureCounts_UCSC:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.star_UCSC.bam",
		ref=config['GTF']['UCSC'],
		script=NGS_PIPELINE + "/scripts/featureCounts.R",
	output:
		gene="{base}/{TIME}/{sample}/TPM_UCSC/{sample}_counts.Gene.txt",
	params:
		rulename   = "featureCounts",
		batch      =config[config['host']]['job_featCount'],
		work_dir =  WORK_DIR
	shell: """
	#######################
	module load R
	cd ${{LOCAL}}
	{input.script} --nt ${{THREADS}} --lib="{wildcards.sample}" --targetFile="{params.work_dir}/{input.bam}" --referenceGTF="{input.ref}" --countOut="{params.work_dir}/{wildcards.base}/{wildcards.TIME}/{wildcards.sample}/TPM_UCSC/{wildcards.sample}_counts" --fpkmOut="{params.work_dir}/{wildcards.base}/{wildcards.TIME}/{wildcards.sample}/TPM_UCSC/{wildcards.sample}_fpkm"
	#######################
	"""

############
## featureCounts
##############
rule FeatureCounts_ENS:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.star_ENS.bam",
		ref=config['GTF']['ENS'],
		script=NGS_PIPELINE + "/scripts/featureCounts.R",
	output:
		"{base}/{TIME}/{sample}/TPM_ENS/{sample}_counts.Exon.fc.RDS",
		"{base}/{TIME}/{sample}/TPM_ENS/{sample}_counts.Exon.txt",
		"{base}/{TIME}/{sample}/TPM_ENS/{sample}_counts.Gene.fc.RDS",
		"{base}/{TIME}/{sample}/TPM_ENS/{sample}_counts.Gene.txt",
		"{base}/{TIME}/{sample}/TPM_ENS/{sample}_counts.Transcript.fc.RDS",
		"{base}/{TIME}/{sample}/TPM_ENS/{sample}_counts.Transcript.txt",
		"{base}/{TIME}/{sample}/TPM_ENS/{sample}_fpkm.Exon.txt",
		"{base}/{TIME}/{sample}/TPM_ENS/{sample}_fpkm.Gene.txt",
		"{base}/{TIME}/{sample}/TPM_ENS/{sample}_fpkm.Transcript.txt"
	params:
		rulename   = "featureCounts",
		batch      =config[config['host']]['job_featCount'],
		work_dir =  WORK_DIR
	shell: """
	#######################
	module load R
	cd ${{LOCAL}}
	{input.script} --nt ${{THREADS}} --lib="{wildcards.sample}" --targetFile="{params.work_dir}/{input.bam}" --referenceGTF="{input.ref}" --countOut="{params.work_dir}/{wildcards.base}/{wildcards.TIME}/{wildcards.sample}/TPM_ENS/{wildcards.sample}_counts" --fpkmOut="{params.work_dir}/{wildcards.base}/{wildcards.TIME}/{wildcards.sample}/TPM_ENS/{wildcards.sample}_fpkm"
	#######################
	"""
