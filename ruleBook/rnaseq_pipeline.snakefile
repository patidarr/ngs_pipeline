RNASEQ_BAM =[]
RNA_QC_ALL =[]
RNASEQ_FUSION =[]
EXPRESSION=[]
RNA_CALLS =[]
SUB2RNA = {}
SUB_RNASEQ=[]
SUB_FUSION={}
SUB_QC={}
for subject,samples in config['RNASeq'].items():
	SUB_RNASEQ.append(subject)
	for sample in samples:
		SUB2RNA[sample]=subject
for subject  in config['RNASeq'].keys():
	DBFiles         +=[subject+"/"+TIME+"/"+subject+"/db/"+subject+".rnaseq"]
	ActionableFiles +=[subject+"/"+TIME+ACT_DIR+subject+".fusion.actionable.txt"]
	ActionableFiles +=[subject+"/"+TIME+ACT_DIR+subject+".rnaseq.actionable.txt"]
	ALL_QC          +=[subject+"/"+TIME+"/qc/"+subject+".RnaSeqQC.txt"]
	RNA_QC_ALL	+=[subject+"/"+TIME+"/qc/"+subject+".RnaSeqQC.txt"]
	ALL_QC		+=[subject+"/"+TIME+"/qc/"+subject+".transcriptCoverage.png"]
	for sample in config['RNASeq'][subject]:
		ALL_QC    +=  [subject+"/"+TIME+"/qc/"+subject+".circos.png"]
		RNASEQ_BAM += [subject+"/"+TIME+"/"+sample+"/"+sample+".star.final.bam"]
		RNASEQ_BAM += [subject+"/"+TIME+"/"+sample+"/"+sample+"_ucsc.SJ.out.tab"]
		RNASEQ_BAM += [subject+"/"+TIME+"/"+sample+"/"+sample+".tophat.final.bam"]
		ALL_FASTQC += [subject+"/"+TIME+"/"+sample+"/qc/fastqc/"+sample+"_R2_fastqc.html"]
		RNASEQ_FUSION += [subject+"/"+TIME+"/"+sample+"/fusion/tophat-fusion.txt"]
		RNASEQ_FUSION += [subject+"/"+TIME+"/"+sample+"/fusion/fusion-catcher.txt"]
		RNASEQ_FUSION += [subject+"/"+TIME+"/"+sample+"/fusion/defuse.filtered.txt"]
		ALL_QC      += [subject+"/"+TIME+"/"+sample+"/qc/"+sample+".star.flagstat.txt"]
		ALL_QC      += [subject+"/"+TIME+"/"+sample+"/qc/"+sample+".star.squeeze.done"]
		ALL_QC      += [subject+"/"+TIME+"/"+sample+"/qc/"+sample+".star.hotspot.depth"]
		ALL_QC      += [subject+"/"+TIME+"/"+sample+"/qc/"+sample+".RnaSeqMetrics.txt"]
		ALL_QC      += [subject+"/"+TIME+"/"+sample+"/qc/"+sample+".RnaSeqMetrics.pdf"]
		ALL_QC      += [subject+"/"+TIME+"/"+sample+"/qc/"+sample+".RnaSeqQC.txt"]
		ALL_QC      += [subject+"/"+TIME+"/"+sample+"/qc/"+sample+".star.gt"]
		ALL_QC      += [subject+"/"+TIME+"/"+sample+"/fusion/"+sample+".actionable.fusion.txt"]
		add_to_SUBJECT_ANNO(subject, "rnaseq", [subject+"/"+TIME+"/"+sample+"/calls/"+sample+".HC_RNASeq.annotated.txt"])
		#EXPRESSION += [subject+"/"+TIME+"/"+sample+"/exonExp_UCSC/"+sample+".exonExpression.UCSC.txt"]
		#EXPRESSION += [subject+"/"+TIME+"/"+sample+"/exonExp_ENS/"+sample+".exonExpression.ENS.txt"]
		for gtf in config['GTF']:
			#EXPRESSION += [subject+"/"+TIME+"/"+sample+"/cufflinks_"+gtf+"/genes.fpkm_tracking_log2"]
			EXPRESSION += [subject+"/"+TIME+"/"+sample+"/TPM_"+gtf+"/"+sample+"_counts.Gene.txt"]
			#ALL_QC     += [subject+"/"+TIME+"/"+subject+"/db/matrixInput_"+sample+"_"+gtf]
	RNA_CALLS  += ["{subject}/{TIME}/{sample}/calls/{sample}.HC_RNASeq.raw.vcf".format(TIME=TIME, subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
	SUB_FUSION[subject] = ["{subject}/{TIME}/{sample}/fusion/{sample}.actionable.fusion.txt".format(TIME=TIME, subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
	SUB_QC[subject]     = ["{subject}/{TIME}/{sample}/qc/{sample}.RnaSeqQC.txt".format(TIME=TIME, subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
	if subject in SUBJECT_VCFS:
		SUBJECT_VCFS[subject] += ["{subject}/{TIME}/{sample}/calls/{sample}.HC_RNASeq.snpEff.txt".format(TIME=TIME, subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
	else:
		SUBJECT_VCFS[subject] = []
		SUBJECT_VCFS[subject] += ["{subject}/{TIME}/{sample}/calls/{sample}.HC_RNASeq.snpEff.txt".format(TIME=TIME, subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
	if subject in SUB_HOT:
		SUB_HOT[subject] += ["{subject}/{TIME}/{sample}/qc/{sample}.star.hotspot.depth".format(TIME=TIME, subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_LOH[subject] += ["{subject}/{TIME}/{sample}/qc/{sample}.star.loh".format(TIME=TIME, subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_COV[subject] += ["{subject}/{TIME}/{sample}/qc/{sample}.star.coverage.txt".format(TIME=TIME, subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_GT[subject]  += ["{subject}/{TIME}/{sample}/qc/{sample}.star.gt".format(TIME=TIME,subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
	else:
		SUB_HOT[subject] = []
		SUB_LOH[subject] = []
		SUB_COV[subject] = []
		SUB_GT[subject]  = []
		SUB_HOT[subject] += ["{subject}/{TIME}/{sample}/qc/{sample}.star.hotspot.depth".format(TIME=TIME, subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_LOH[subject] += ["{subject}/{TIME}/{sample}/qc/{sample}.star.loh".format(TIME=TIME, subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_COV[subject] += ["{subject}/{TIME}/{sample}/qc/{sample}.star.coverage.txt".format(TIME=TIME,subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_GT[subject]  += ["{subject}/{TIME}/{sample}/qc/{sample}.star.gt".format(TIME=TIME, subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
############
#       RNASeq All
############
rule RNASeq:
	input:
		RNASEQ_BAM,
		EXPRESSION,
		RNA_CALLS,
		RNASEQ_FUSION,
		expand("{subject}/{TIME}"+ACT_DIR+"{subject}.fusion.actionable.txt", TIME=TIME, subject=config['RNASeq']),
	output:
		expand("ngs_pipeline_{NOW}.rnaseq.done", NOW=NOW)
	wildcard_constraints:
                NOW="\w+"	
	params:
		rulename  = "RNASeq_final",
		batch     = config[config['host']]["job_default"]

	shell: """
	#######################
	touch {output}
	#######################
	"""
############
#	Tophat
############
rule TOPHAT:
	input: R=lambda wildcards: FQ[wildcards.sample],
	output:
		"{base}/{TIME}/{sample}/tophat_{sample}/accepted_hits.bam",
		"{base}/{TIME}/{sample}/tophat_{sample}/accepted_hits.bam.bai"
	version: config["tophat"]
	params:
		rulename  = "tophat",
		samtools  = config['samtools'],
		batch     = config[config['host']]["job_tophat"],
		ref=config['Bowtie2Index']
	shell: """
	#######################
	module load tophat/{version}
	module load samtools/{params.samtools}
	tophat -p ${{THREADS}} -o ${{LOCAL}} --keep-fasta-order --rg-id {wildcards.sample} --no-coverage-search --rg-sample {wildcards.sample} --rg-library {wildcards.sample} --rg-platform ILLUMINA --fusion-search --fusion-min-dist 100000 --mate-inner-dist 84 --mate-std-dev 74 {params.ref} {input.R[0]} {input.R[1]}
	cp -rf ${{LOCAL}}/* {wildcards.base}/{TIME}/{wildcards.sample}/tophat_{wildcards.sample}/
	samtools index {wildcards.base}/{TIME}/{wildcards.sample}/tophat_{wildcards.sample}/accepted_hits.bam
	#######################
	"""
############
#       Link Tophat bam file
############
rule TOPHAT_LINK:
	input:
		bam="{base}/{TIME}/{sample}/tophat_{sample}/accepted_hits.bam",
		bai="{base}/{TIME}/{sample}/tophat_{sample}/accepted_hits.bam.bai",
	output:
		bam="{base}/{TIME}/{sample}/{sample}.tophat.final.bam",
		bai="{base}/{TIME}/{sample}/{sample}.tophat.final.bam.bai"
	params:
		rulename  = "tophat_link",
		batch     = config[config['host']]["job_default"]
	shell: """
	#######################
	sleep 60
	cd {wildcards.base}/{TIME}/{wildcards.sample}/
	ln -sf tophat_{wildcards.sample}/accepted_hits.bam {wildcards.sample}.tophat.final.bam
	ln -sf tophat_{wildcards.sample}/accepted_hits.bam.bai {wildcards.sample}.tophat.final.bam.bai
	#######################
	"""
############
#       Tophat-fusion
############
rule TOPHAT_FUSION:
	input:
		"{base}/{TIME}/{sample}/tophat_{sample}/accepted_hits.bam",
		"{base}/{TIME}/{sample}/tophat_{sample}/accepted_hits.bam.bai"
	output:
		"{base}/{TIME}/{sample}/tophatfusion_out/result.txt",
		"{base}/{TIME}/{sample}/fusion/tophat-fusion.txt"
	version: config["tophat"]
	params:
		rulename = "tp",
		blast	 =config['version_blast'],
		batch    =config[config['host']]['job_tophatPost'],
		ref      =config['BowtieIndex'],
		bowtie   =config['bowtie'],
		tp_ref   =config['tophat_post_ref']
	shell: """
	#######################
	module load tophat/{version}
	module load bowtie/{params.bowtie}
	module load blast/{params.blast}
	cd {wildcards.base}/{TIME}/{wildcards.sample}/
	rm -f blast ensGene.txt ensGtp.txt mcl refGene.txt
	ln -s {params.tp_ref}/* .
	tophat-fusion-post -p ${{THREADS}} --num-fusion-pairs 1 {params.ref}
	rm blast ensGene.txt ensGtp.txt mcl refGene.txt
	sed -i  '1s/^/Sample\\tGene_left\\tChr_left\\tCoordinate_left\\tGene_right\\tChr_right\\tCoordinate_right\\t#SpanningReads\\t#SpanningMatePairs\\t#SpanningMateEndOfPair\\tScore\\n/' tophatfusion_out/result.txt
	ln -sf ../tophatfusion_out/result.html fusion/tophat-fusion.html
	ln -sf ../tophatfusion_out/result.txt  fusion/tophat-fusion.txt
	#######################
	"""
############
#       Cufflinks
############
rule CUFFLINKS:
	input:
		bam="{base}/{TIME}/{sample}/tophat_{sample}/accepted_hits.bam",
		bai="{base}/{TIME}/{sample}/tophat_{sample}/accepted_hits.bam.bai",
		convertor =NGS_PIPELINE + "/scripts/fpkm2log2_fpkm.pl",
		ref=lambda wildcards: config['GTF'][wildcards.gtf]
	output:
		gene="{base}/{TIME}/{sample}/cufflinks_{gtf}/genes.fpkm_tracking",
		gene_log="{base}/{TIME}/{sample}/cufflinks_{gtf}/genes.fpkm_tracking_log2",
		isoform="{base}/{TIME}/{sample}/cufflinks_{gtf}/isoforms.fpkm_tracking",
		isoform_log="{base}/{TIME}/{sample}/cufflinks_{gtf}/isoforms.fpkm_tracking_log2"
	version: config['cufflinks']
	params:
		rulename   = "cuff",
		batch      =config[config['host']]['job_cufflinks']
	shell: """
	#######################
	module load cufflinks/{version}
	cufflinks --no-update-check -p ${{THREADS}} -G {input.ref} --max-bundle-frags 8000000000000 --max-bundle-length 10000000 -o {wildcards.base}/{TIME}/{wildcards.sample}/cufflinks_{wildcards.gtf} {input.bam}
	perl {input.convertor} {output.gene}    > {output.gene_log}
	perl {input.convertor} {output.isoform} > {output.isoform_log}
	#######################
	"""
############
#	This is to make input file for Sivasish's metrix generation tool.
############
rule Cuff_Mat:
	input: 
		gene="{base}/{TIME}/{sample}/cufflinks_{gtf}/genes.fpkm_tracking"
	output:
		"{base}/{TIME}/{base}/db/matrixInput_{sample}_{gtf}"
	params:
		rulename  ="cuff_1",
		batch     =config[config['host']]['job_default'],
		diagnosis =lambda wildcards: config['Diagnosis'][wildcards.sample],
		library   =lambda wildcards: config['sample_captures'][wildcards.sample]
	shell: """
	#######################
	echo -e "{wildcards.sample}\\t{params.diagnosis}\\t{params.library}" >{output}
	echo "{WORK_DIR}/{input.gene}" >>{output}
	#######################
	"""
############
#       Exon Expression
############
rule EXON_EXP_UCSC:
	input:
		bam="{base}/{TIME}/{sample}/tophat_{sample}/accepted_hits.bam",
		bai="{base}/{TIME}/{sample}/tophat_{sample}/accepted_hits.bam.bai",
		convertor =NGS_PIPELINE + "/scripts/exon_exp.sh",
		ucsc=config["exon_Bed_UCSC"]
	output:
		ucsc="{base}/{TIME}/{sample}/exonExp_UCSC/{sample}.exonExpression.UCSC.txt"
	version: config['samtools']
	params:
		rulename   = "exonExp.1",
		batch      =config[config['host']]['job_exonExp']
	shell: """
	#######################
	module load samtools/{version}
	totalReads=`samtools flagstat {input.bam} |head -1 | sed 's/\s/\\t/g' | cut -f1`

        split -d -l 20000 {input.ucsc} ${{LOCAL}}/ucsc
        for file in ${{LOCAL}}/ucsc*
        do
		sh {input.convertor} ${{totalReads}} ${{file}} {input.bam} ${{file}}.out &
        done
        wait;
        cat ${{LOCAL}}/ucsc*.out >{output.ucsc}
	rm -rf ${{LOCAL}}/*
	#######################
        """
############
#       Exon Expression
############
rule EXON_EXP_ENS:
	input:
		bam="{base}/{TIME}/{sample}/tophat_{sample}/accepted_hits.bam",
		bai="{base}/{TIME}/{sample}/tophat_{sample}/accepted_hits.bam.bai",
		convertor =NGS_PIPELINE + "/scripts/exon_exp.sh",
		ens=config["exon_Bed_ENS"],
	output:
		ens="{base}/{TIME}/{sample}/exonExp_ENS/{sample}.exonExpression.ENS.txt"
	version: config['samtools']
	params:
		rulename   = "exonExp.2",
		batch      =config[config['host']]['job_exonExp']
	shell: """
	#######################
	module load samtools/{version}
	totalReads=`samtools flagstat {input.bam} |head -1 | sed 's/\s/\\t/g' | cut -f1`
	
	split -d -l 30000 {input.ens} ${{LOCAL}}/ens
	for file in ${{LOCAL}}/ens*
	do
		sh {input.convertor} ${{totalReads}} ${{file}} {input.bam} ${{file}}.out &
	done
	wait
	cat ${{LOCAL}}/ens*.out >{output.ens}
	#######################
	"""
############
#       Fusioncatcher
############
rule FUSION_CATCHER:
	input: R=lambda wildcards: FQ[wildcards.sample]
	output:
		"{base}/{TIME}/{sample}/fusion/fusion-catcher.txt"
	version: config['fusioncatcher']
	params:
		rulename = "fc",
		batch    = config[config['host']]['job_fusioncatch']
	shell: """
	#######################
	module load fusioncatcher/{version}
	fusioncatcher -p ${{THREADS}} -i {input.R[0]},{input.R[1]} -o ${{LOCAL}}/
	cp ${{LOCAL}}/final-list_candidate-fusion-genes.GRCh37.txt {wildcards.base}/{TIME}/{wildcards.sample}/fusion/fusion-catcher.txt
	#######################
	"""
############
#       deFuse
############
rule DeFuse:
	input: R=lambda wildcards: FQ[wildcards.sample],
	output:
		"{base}/{TIME}/{sample}/fusion/defuse.raw.txt",
		"{base}/{TIME}/{sample}/fusion/defuse.filtered.txt",
		"{base}/{TIME}/{sample}/fusion/defuse.Reads/defuse.done"
	version: config["defuse"]
	resources: DeFuse=1
	params:
		rulename = "deFuse",
		batch    = config[config['host']]["job_deFuse"],
		defuse_config=NGS_PIPELINE + "/Tools_config/"+config["defuse_config"],
		parallel = NGS_PIPELINE + "/scripts/parallel"
	shell: """
	#######################
	module load defuse/{version}
	export TMPDIR="${{LOCAL}}/temp/R"
	export TMP="${{LOCAL}}/temp/R"
	export TEMP="${{LOCAL}}/temp/R"

	gunzip -c {input.R[0]} >${{LOCAL}}/{wildcards.sample}_R1.fastq &
	gunzip -c {input.R[1]} >${{LOCAL}}/{wildcards.sample}_R2.fastq &
	wait
	defuse.pl -c {params.defuse_config} \
		-1 ${{LOCAL}}/{wildcards.sample}_R1.fastq\
		-2 ${{LOCAL}}/{wildcards.sample}_R2.fastq\
		-p ${{THREADS}} \
		-n {wildcards.sample}\
		-o ${{LOCAL}}/defuse\
		-s direct
	cp ${{LOCAL}}/defuse/results.filtered.tsv  {wildcards.base}/{TIME}/{wildcards.sample}/fusion/defuse.filtered.txt
	cp ${{LOCAL}}/defuse/results.tsv           {wildcards.base}/{TIME}/{wildcards.sample}/fusion/defuse.raw.txt
	
	for ID in `cut -f1 {wildcards.base}/{TIME}/{wildcards.sample}/fusion/defuse.filtered.txt|grep -v cluster_id`;
	do
		echo "get_reads.pl -c {params.defuse_config} -o ${{LOCAL}}/defuse/ -i ${{ID}} >{wildcards.base}/{TIME}/{wildcards.sample}/fusion/defuse.Reads/${{ID}}.txt"
	done >{wildcards.base}/{TIME}/{wildcards.sample}/fusion/cmd.swarm
	cat {wildcards.base}/{TIME}/{wildcards.sample}/fusion/cmd.swarm | {params.parallel} -j ${{THREADS}} --no-notice
	touch {wildcards.base}/{TIME}/{wildcards.sample}/fusion/defuse.Reads/defuse.done
	rm -rf ${{LOCAL}}/defuse {wildcards.base}/{TIME}/{wildcards.sample}/fusion/cmd.swarm
	#######################
	"""
############
#       STAR
############
rule STAR:
	input:  R=lambda wildcards: FQ[wildcards.sample],
		ref=config["reference"],
	output:
		temp("{base}/{TIME}/{sample}/{sample}.star.bam"),
		temp("{base}/{TIME}/{sample}/{sample}.star.bam.bai")
	version: config["STAR"]
	params:
		rulename  = "STAR",
		samtools  = config['samtools'],
		batch     = config[config['host']]['job_STAR'],
		star_ref  = config['STAR_ref'],
		awk       = NGS_PIPELINE + "/scripts/SJDB.awk",
		home      = WORK_DIR,
		picard    = config['picard']
	shell: """
	#######################
	module load STAR/{version}
	cd ${{LOCAL}}/
	# run 1st pass
	STAR --outTmpDir STEP1 \
		--genomeDir {params.star_ref} \
		--readFilesIn {input.R[0]} {input.R[1]} \
		--readFilesCommand zcat\
		--outFileNamePrefix {wildcards.sample} \
		--runThreadN ${{THREADS}} \
		--outFilterMismatchNmax 2
	echo "Finished Step 1"

	# make splice junctions database file out of SJ.out.tab, filter out non-canonical junctions
	mkdir GenomeForPass2
	awk -f {params.awk} {wildcards.sample}SJ.out.tab > GenomeForPass2/{wildcards.sample}.out.tab.Pass1.sjdb
	echo "Finished Step 2"

	# generate genome with junctions from the 1st pass
	STAR --outTmpDir STEP1\
		--genomeDir GenomeForPass2\
		--runMode genomeGenerate\
		--genomeSAindexNbases 8\
		--genomeFastaFiles {input.ref}\
		--sjdbFileChrStartEnd GenomeForPass2/{wildcards.sample}.out.tab.Pass1.sjdb\
		--sjdbOverhang 100\
		--runThreadN ${{THREADS}}
	echo "Finished Step 3"

	# run 2nd pass with the new genome
	STAR --outTmpDir STEP1\
		--genomeDir GenomeForPass2\
		--runThreadN ${{THREADS}}\
		--outSAMattributes All\
		--readFilesIn {input.R[0]} {input.R[1]}\
		--readFilesCommand zcat\
		--genomeLoad NoSharedMemory\
		--outFileNamePrefix {wildcards.sample}_pass2
	echo "Finished Step 4"

	module load picard/{params.picard}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $PICARD_JAR AddOrReplaceReadGroups\
	VALIDATION_STRINGENCY=SILENT\
	INPUT={wildcards.sample}_pass2Aligned.out.sam\
	OUTPUT={params.home}/{wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}.star.bam\
	SORT_ORDER=coordinate RGLB={wildcards.sample} RGPU={wildcards.sample} RGPL=ILLUMINA RGSM={wildcards.sample} RGCN=khanlab

	echo "Finished Step 5"
	module load samtools/{params.samtools}
	samtools index {params.home}/{wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}.star.bam
	#######################
	"""
############
# RNASeq Hapcaller
############
rule HapCall_RNASeq:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.star.final.bam",
		bai="{base}/{TIME}/{sample}/{sample}.star.final.bam.bai",
		ref=config["reference"],
		dbsnp=config["dbsnp"]
	output:
		vcf="{base}/{TIME}/{sample}/calls/{sample}.HC_RNASeq.raw.vcf"
	version: config["GATK"]
	params:
		rulename = "HC_RNA",
		batch    = config[config['host']]["job_HC"]
	shell: """
	#######################
	module load GATK/{version}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T HaplotypeCaller -R {input.ref} -I {input.bam} -o ${{LOCAL}}/{wildcards.sample}.vcf --dbsnp {input.dbsnp} -dontUseSoftClippedBases -stand_call_conf 30 -stand_emit_conf 30
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T VariantFiltration -R {input.ref} -V ${{LOCAL}}/{wildcards.sample}.vcf -window 35 -cluster 3 --filterExpression "FS > 30.0 || QD < 2" -filterName "RNASeqFilters_FS_QD" -o {output.vcf}
	#######################
	"""
############
# Coverage
############
rule CoveragE:
	input:
		bam="{subject}/{TIME}/{sample}/{sample}.star.final.bam",
		bai="{subject}/{TIME}/{sample}/{sample}.star.final.bam.bai",
		interval=lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]]
	output:
		"{subject}/{TIME}/{sample}/qc/{sample}.star.coverage.txt"
	version: config["bedtools"]
	params:
		rulename = "coverage",
		batch    = config[config['host']]["job_bedtools"]
	shell: """
	#######################
	module load bedtools/{version}
	bedtools coverage -abam {input.bam} -b {input.interval} -hist |grep "^all" > {output}
	#######################
	"""
############
# Filter fusion for every library
############
rule Sub_Fusion:
	input:
		tophat="{subject}/{TIME}/{sample}/fusion/tophat-fusion.txt",
		fc="{subject}/{TIME}/{sample}/fusion/fusion-catcher.txt",
		defuse="{subject}/{TIME}/{sample}/fusion/defuse.filtered.txt",
		convertor = NGS_PIPELINE + "/scripts/" + config['Actionable_fusion'],
	output:
		"{subject}/{TIME}/{sample}/fusion/{sample}.actionable.fusion.txt"
	params:
		rulename = "Sub_F",
		batch    = config[config['host']]["job_default"]
	shell: """
	#######################
	mkdir -p {wildcards.subject}/{TIME}/Actionable
	perl {input.convertor} {wildcards.sample} {input.defuse} {input.tophat} {input.fc} {wildcards.subject}/{TIME}{ACT_DIR} |awk 'NR<2{{print $0;next}}{{print $0| "sort "}}' >{output}
	#######################
	"""
############
# Combine filtered fusions to actionable.
############
rule Actionable_fusion:
	input:
		fusion=lambda wildcards: SUB_FUSION[wildcards.subject]
	output:
		"{subject}/{TIME}/{ACT_DIR}{subject}.fusion.actionable.txt"
	params:
		rulename = "Actionable_F",
		batch    = config[config['host']]["job_default"]
	shell: """
	#######################
	cat {input.fusion} |sort |uniq >{output}.tmp
	grep "LeftGene" {output}.tmp >{output}
	grep -v "LeftGene" {output}.tmp >>{output}
	
	rm -rf {output}.tmp
	#######################
	"""
