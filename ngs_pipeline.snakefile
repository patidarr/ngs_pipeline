import itertools
import os
import collections
import json
from snakemake.utils import R
from snakemake.exceptions import MissingInputException
# Snakemake Base location

NGS_PIPELINE=os.environ['NGS_PIPELINE']
WORK_DIR=os.environ['WORK_DIR']
DATA_DIR=os.environ['DATA_DIR']
ACT_DIR=os.environ['ACT_DIR']
HOST=os.environ['HOST']
TIME=os.environ['TIME']
if HOST == 'biowulf.nih.gov':
	configfile: NGS_PIPELINE +"/config/config_common.json"
	configfile: NGS_PIPELINE +"/config/config_common_biowulf.json"
	configfile: NGS_PIPELINE +"/config/config_cluster.json"
elif HOST == 'login01':
	configfile: NGS_PIPELINE +"/config/config_common.json"
	configfile: NGS_PIPELINE +"/config/config_common_tgen.json"
	configfile: NGS_PIPELINE +"/config/config_cluster.json"

config['host'] = HOST
#HOST = config['host']
###########################################################################
#
#		This initializes all the variables we need for the jobs.
#		It also removes the host specific constraints like scratch
#		area on the node.
#		module purge is needed to remove all the loaded modules and
#			inside the rule load what is necessary.
###########################################################################
shell.prefix("""
set -e -o pipefail
module purge
if [ {HOST} == 'biowulf.nih.gov' ]
	then
		MEM=`echo "${{SLURM_MEM_PER_NODE}} / 1024 "|bc`
		LOCAL="/lscratch/${{SLURM_JOBID}}/"
		THREADS=${{SLURM_CPUS_ON_NODE}}
elif [ {HOST} == 'login01' ]
	then
		module load torque/4.2.2
		module load gcc/4.8.1
		MEM=`qstat -f ${{PBS_JOBID}} |grep Resource_List.mem |perl -n -e'/mem = (\d+)gb/ && print \$1'`
		mkdir -p /projects/scratch/${{PBS_JOBID}}/
		LOCAL="/projects/scratch/${{PBS_JOBID}}/"
		THREADS=${{PBS_NUM_PPN}}
fi
""")
###########################################################################
#
#			Conversion
#
###########################################################################
SUBJECT_TO_SAMPLE  = {}
for subject in config['subject']:
	SUBJECT_TO_SAMPLE[subject] = expand("{sample}", sample = config['subject'][subject])
###########################################################################
SAMPLE_TO_SUBJECT  = {}
for subject,samples in config['subject'].items():
	for sample in samples:
		SAMPLE_TO_SUBJECT[sample]=subject
###########################################################################
####
#### Targets
####
PATIENTS =[]
SUBS  = []
SUB_BAMS= {}
SUB_COV = {}
SUB_LOH = {}
SUB_GT  = {}
SUB_MPG = {}
SUB_HOT = {}
SUB_IGV = {}
SUB_CON_QC = {}
SAMPLES =[]
somaticPairs = {}
somaticCopy = {}
pairedCapture = {}
FQ={}

#for sub in config["subject"]:
for sample in config['library'].keys():
	for fq in config['library'][sample]:
		if len(config['library'][sample]) == 1:
			FQ[sample] =[DATA_DIR+fq+"/"+fq+"_R1.fastq.gz", DATA_DIR+fq+"/"+fq+"_R2.fastq.gz"]
		else:
			exit()

for subject in config['subject'].keys():
	SUBS.append(subject)
	PATIENTS.append(subject)
	SUB_BAMS[subject]= ["{subject}/{TIME}/{sample}/{sample}.bwa.final.bam".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_COV[subject] = ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.coverage.txt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_HOT[subject] = ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.hotspot.depth".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_LOH[subject] = ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.loh".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_GT[subject]  = ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.gt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_CON_QC[subject]  = ["{subject}/{TIME}/{sample}/qc/{sample}.consolidated_QC".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_MPG[subject] = ["{subject}/{TIME}/{sample}/calls/{sample}.bam2mpg.vcf.gz".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_IGV[subject] = ["{subject}/{TIME}/{sample}/{sample}.bwa.final.bam".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_IGV[subject]+= ["{subject}/{TIME}/{sample}/{sample}.bwa.final.bam.tdf".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_IGV[subject]+= ["{subject}/{TIME}/{sample}/{sample}.novo.final.bam".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_IGV[subject]+= ["{subject}/{TIME}/{sample}/{sample}.novo.final.bam.tdf".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	for sample in config['subject'][subject]:
		SAMPLES.append(sample)
###########################################################################
#		Add RNASeq only samples to PATIENTS
###########################################################################
for subject in config['RNASeq']:
	SUBJECT_TO_SAMPLE[subject] = expand("{sample}", sample = config['RNASeq'][subject])
for subject  in config['RNASeq'].keys():
        if subject not in PATIENTS:
                PATIENTS.append(subject)
###########################################################################
ALL_FASTQC  = ["{subject}/{TIME}/{sample}/qc/fastqc/{sample}_R2_fastqc.html".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC      = ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.flagstat.txt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC     += ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.hotspot.depth".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC     += ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.gt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC     += ["{subject}/{TIME}/{sample}/qc/BamQC/qualimapReport.html".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC	   += ["{subject}/{TIME}/{sample}/qc/{sample}.depth_per_base".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC	   += ["{subject}/{TIME}/{sample}/qc/{sample}.failExons".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC	   += ["{subject}/{TIME}/{sample}/qc/{sample}.failGenes".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC	   += ["{subject}/{TIME}/{sample}/qc/{sample}.hsmetrics".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
CON_QC	    = ["{subject}/{TIME}/{sample}/qc/{sample}.consolidated_QC".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC     += ["{subject}/{TIME}/{sample}/copyNumber/{sample}.count.txt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC     += expand("{subject}/{TIME}/qc/{subject}.genotyping.txt", TIME=TIME, subject=PATIENTS)
ALL_QC     += expand("{subject}/{TIME}/qc/{subject}.consolidated_QC.txt", TIME=TIME, subject=SUBS)
ALL_QC     += expand("{subject}/{TIME}/qc/{subject}.coveragePlot.png",TIME=TIME, subject=PATIENTS)
ALL_QC     += expand("{subject}/{TIME}/qc/{subject}.circos.png", TIME=TIME, subject=PATIENTS)
ALL_QC     += expand("{subject}/{TIME}/qc/{subject}.hotspot_coverage.png", TIME=TIME, subject=PATIENTS)
ALL_QC     += expand("{subject}/annotation/{subject}.Annotations.coding.rare.txt", subject=PATIENTS)
ALL_QC     += expand("{subject}/{TIME}/igv/session_{subject}.xml", TIME=TIME, subject=PATIENTS)
if len(config['sample_references']) > 0:
	for Tumor in config['sample_references']:
		for Normal in config['sample_references'][Tumor]:
			TumorBam   = "{subject}/{TIME}/{sample}/{sample}.bwa.final".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[Tumor], sample=Tumor)
			TumorCopy  = "{subject}/{TIME}/{sample}/copyNumber/{sample}.count.txt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[Tumor], sample=Tumor)
			NormalBam  = "{subject}/{TIME}/{sample}/{sample}.bwa.final".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[Normal], sample=Normal)
			NormalCopy = "{subject}/{TIME}/{sample}/copyNumber/{sample}.count.txt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[Normal], sample=Normal)
			pairedCapture[Tumor] = config['sample_captures'][Tumor]
			somaticPairs[Tumor] = [TumorBam + ".bam" , TumorBam + ".bam.bai", NormalBam + ".bam", NormalBam + ".bam.bai"]
			somaticCopy[Tumor] = [NormalCopy, TumorCopy]
###########################################################################
SUBJECT_ANNO = dict([(key, {}) for key in PATIENTS])
def add_to_SUBJECT_ANNO(subject, category, file_list):
	if category not in SUBJECT_ANNO[subject]:
		SUBJECT_ANNO[subject][category] = file_list
	else:
		SUBJECT_ANNO[subject][category].extend(file_list)
###########################################################################
SUBJECT_VCFS = {}
COPY_NUMBER=[]
SOMATIC =[]
#for subject in SUBS:
#	local  = []
#	local.extend([(subject+"/"+subject+"/calls/"+subject+".hapCaller.snpEff.txt"),
#		      (subject+"/"+subject+"/calls/"+subject+".platypus.snpEff.txt"),
#		      (subject+"/"+subject+"/calls/"+subject+".bam2mpg.snpEff.txt")])
#	if subject not in SUBJECT_VCFS:
#		SUBJECT_VCFS[subject] = local
#	germline = [w.replace('snpEff','annotated') for w in local]
#	add_to_SUBJECT_ANNO(subject,"germline",germline)

for subject in config['subject']:
	local  = []
	for sample in config['subject'][subject]:
		local.extend([(subject+"/"+TIME+"/"+sample+"/calls/"+sample+".HC_DNASeq.snpEff.txt"),
			      (subject+"/"+TIME+"/"+sample+"/calls/"+sample+".Platypus.snpEff.txt"),
			      (subject+"/"+TIME+"/"+sample+"/calls/"+sample+".bam2mpg.snpEff.txt")])
	if subject not in SUBJECT_VCFS:
		SUBJECT_VCFS[subject] = local
	germline = [w.replace('snpEff','annotated') for w in local]
	add_to_SUBJECT_ANNO(subject,"germline",germline)
	

for sample in config['sample_references'].keys():
	local  = []
	subject=SAMPLE_TO_SUBJECT[sample]
	local.extend(
		[ (subject+"/"+TIME+"/"+sample+"/calls/"+sample+".MuTect.snpEff.txt"),
		  (subject+"/"+TIME+"/"+sample+"/calls/"+sample+".strelka.snvs.snpEff.txt"),
		  (subject+"/"+TIME+"/"+sample+"/calls/"+sample+".strelka.indels.snpEff.txt")
		]
	)
	COPY_NUMBER +=[subject+"/"+TIME+"/"+sample+"/copyNumber/"+sample+".copyNumber.txt"]
	COPY_NUMBER +=[subject+"/"+TIME+"/"+sample+"/copyNumber/"+sample+".hq.txt"]
	COPY_NUMBER +=[subject+"/"+TIME+"/"+sample+"/copyNumber/"+sample+".CN.annotated.txt"]
	COPY_NUMBER +=[subject+"/"+TIME+"/"+sample+"/copyNumber/"+sample+".CN.filtered.txt"]
	COPY_NUMBER +=[subject+"/"+TIME+"/"+sample+"/copyNumber/"+sample+".CN.png"]
	SOMATIC     +=[subject+"/"+TIME+"/"+sample+"/calls/"+sample+".MuTect.annotated.txt"]
	SOMATIC     +=[subject+"/"+TIME+"/"+sample+"/calls/"+sample+".strelka.snvs.annotated.txt"]
	SOMATIC     +=[subject+"/"+TIME+"/"+sample+"/calls/"+sample+".strelka.indels.annotated.txt"]
	if subject in SUBJECT_VCFS:
		SUBJECT_VCFS[subject].extend(local)

	somatic = [w.replace('snpEff','annotated') for w in local]
	if sample in config['sample_RNASeq']:
		somatic = [w.replace('MuTect.annotated','MuTect.annotated.expressed') for w in somatic]
		somatic = [w.replace('strelka.snvs.annotated','strelka.snvs.annotated.expressed') for w in somatic]
		somatic = [w.replace('strelka.indels.annotated','strelka.indels.annotated.expressed') for w in somatic]
	add_to_SUBJECT_ANNO(subject,"somatic",somatic)
###########################################################################
###########################################################################
ALL_EXPRESSED =[]
expressedPairs = {}
if len(config['sample_RNASeq']) > 0:
	for Tumor in config['sample_RNASeq']:
		for RNASample in config['sample_RNASeq'][Tumor]:
			subject=SAMPLE_TO_SUBJECT[Tumor]
			RNASeqBam    = subject + "/"+TIME+ "/" + RNASample + "/calls/"+RNASample + ".HC_RNASeq.snpEff.txt"
			expressedPairs[Tumor] = RNASeqBam
			ALL_EXPRESSED += ["{subject}/{TIME}/{sample}/calls/{sample}.MuTect.annotated.expressed.txt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[Tumor],  sample=Tumor)]
			ALL_EXPRESSED += ["{subject}/{TIME}/{sample}/calls/{sample}.strelka.snvs.annotated.expressed.txt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[Tumor], sample=Tumor)]
			ALL_EXPRESSED += ["{subject}/{TIME}/{sample}/calls/{sample}.strelka.indels.annotated.expressed.txt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[Tumor], sample=Tumor)]
###########################################################################
# we have to do it this way as some samples may not have rna or tumor     #
###########################################################################
varFiles=[]
DBFiles =[]
ActionableFiles =[]
for subject in SUBJECT_ANNO.keys():
	for group in SUBJECT_ANNO[subject].keys():
		DBFiles +=[subject+"/"+TIME+"/"+subject+"/db/"+subject+"."+group]
		ActionableFiles +=[subject+"/"+TIME+ACT_DIR+subject+"."+group+".actionable.txt"]
		for varFile in SUBJECT_ANNO[subject][group]:
			varFiles.append(varFile)
###########################################################################
localrules: Khanlab_Pipeline,
#RNASeq, IGV_Session, DBinput, AttachAnnotation, Expressed, vcf2txt, symlink_tophatBam, copyNovoBam, Actionable_Germline, Actionable_RNAseq, Actionable_Somatic, Circos, CoveragePlot, BoxPlot_Hotspot
###########################################################################
#                               RNASeq Rules                              #
###########################################################################
include: NGS_PIPELINE +"/ruleBook/rnaseq_pipeline.snakefile"
include: NGS_PIPELINE +"/ruleBook/readDepth.snakefile"
include: NGS_PIPELINE +"/ruleBook/failedExon.snakefile"
include: NGS_PIPELINE +"/ruleBook/hsMetrix.snakefile"
include: NGS_PIPELINE +"/ruleBook/Consolidate.snakefile"
include: NGS_PIPELINE +"/ruleBook/universal.snakefile"
include: NGS_PIPELINE +"/ruleBook/haplotypeCaller.snakefile"
include: NGS_PIPELINE +"/ruleBook/platypus.snakefile"
include: NGS_PIPELINE +"/ruleBook/gatk_RNASeq.snakefile"
include: NGS_PIPELINE +"/ruleBook/ideogram.snakefile"
ALL_VCFs =[]
for subject in SUBJECT_VCFS.keys():
	for vcf in SUBJECT_VCFS[subject]:
		vcf = vcf.replace('snpEff.txt', 'raw.vcf')
		ALL_VCFs +=[vcf]
		vcf = vcf.replace('raw.vcf', 'raw.snpEff.vcf')
		ALL_VCFs +=[vcf]
rule Khanlab_Pipeline:
	input:
		SUB_IGV.values(),
		COPY_NUMBER,
		ALL_VCFs,
		"rnaseqDone",
		CON_QC,
		"Consolidated_QC.txt",
		ALL_QC,
		ALL_FASTQC,
		varFiles,
		DBFiles,
		ActionableFiles
	version: "1.0"
	params:
		rulename = "Final",
		group    = config["group"],
		wait4job = NGS_PIPELINE + "/scripts/block_for_jobid.pl",
		sort 	 = NGS_PIPELINE + "/scripts/awk_sort_withHeader.awk",
		mail 	 = NGS_PIPELINE + "/scripts/tsv2html.final.sh",
		email     = config["mail"],
		host     = config["host"],
		subs     = config["subject"].keys()
	shell: """
	#######################
	find log/ -type f -empty -delete
	find . -group $USER -exec chgrp {params.group} {{}} \;
	find . -type f -user $USER -exec chmod g+r {{}} \; 
#	find . \( -type f -user $USER -exec chmod g+r {{}} \; \) , \( -type d -user $USER -exec chmod g+rwxs {{}} \; \)
	cut -f 1,3 {WORK_DIR}/Consolidated_QC.txt |{params.sort} |uniq >{WORK_DIR}/tmpFile.txt
	ssh {params.host} "{params.mail} --location {WORK_DIR} --host {params.host} --head {WORK_DIR}/tmpFile.txt |/usr/bin/mutt -e \\\"my_hdr Content-Type: text/html\\\" -s 'Khanlab ngs-pipeline Status' `whoami`@mail.nih.gov {params.email}"
	rm -rf {WORK_DIR}/tmpFile.txt
	#######################
	"""
############
#	FASTQC
############
rule FASTQC:
	input: R=lambda wildcards: FQ[wildcards.sample]
	output:
		"{base}/{TIME}/{sample}/qc/fastqc/{sample}_R1_fastqc.html",
		"{base}/{TIME}/{sample}/qc/fastqc/{sample}_R2_fastqc.html"
	version: config["fastqc"]
	params:
		rulename  = "fastqc",
		batch     = config[config['host']]["job_fastqc"]
	shell: """
	#######################
	module load fastqc/{version}
	ln -sf {input.R[0]} ${{LOCAL}}/{wildcards.sample}_R1.fastq.gz
	ln -sf {input.R[1]} ${{LOCAL}}/{wildcards.sample}_R2.fastq.gz

	fastqc --extract -t ${{THREADS}} -o {wildcards.base}/{TIME}/{wildcards.sample}/qc/fastqc/ -d ${{LOCAL}} ${{LOCAL}}/{wildcards.sample}_R1.fastq.gz
	fastqc --extract -t ${{THREADS}} -o {wildcards.base}/{TIME}/{wildcards.sample}/qc/fastqc/ -d ${{LOCAL}} ${{LOCAL}}/{wildcards.sample}_R2.fastq.gz
	#######################
	"""
############
#       BWA
############
rule BWA:
	input: R=lambda wildcards: FQ[wildcards.sample],
		ref=config["bwaIndex"]
	output:
		temp("{base}/{TIME}/{sample}/{sample}.bwa.bam"),
		temp("{base}/{TIME}/{sample}/{sample}.bwa.bam.bai")
	version: config["bwa"]
	params:
		rulename  = "bwa",
		platform  = config["platform"],
		samtools  = config["samtools"],
		batch     = config[config['host']]["job_bwa"]
	shell: """
	#######################
	module load bwa/{version}
	module load samtools/{params.samtools}
	bwa mem -M \
	-t ${{THREADS}}\
	-R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:{params.platform}' \
	{input.ref} {input.R[0]} {input.R[1]} | samtools view -Sbh - \
	| samtools sort -m 30000000000 - {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}.bwa
	samtools index {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}.bwa.bam
	#######################
	"""
############
#	Novoalign
############
rule NOVOALIGN:
	input: R=lambda wildcards: FQ[wildcards.sample],
		index=config["Novo_index"]
	output:
		temp("{subject}/{TIME}/{sample}/{sample}.novo.bam"),
		temp("{subject}/{TIME}/{sample}/{sample}.novo.bam.bai")
	version: config["novocraft"]
	resources: novoalign=1
	params:
		rulename  = "novoalign",
		batch     = config[config['host']]["job_novoalign"],
		samtools  = config["samtools"],
		platform  = config["platform"]
	shell: """
	#######################
	module load samtools/{params.samtools}
	module load novocraft/{version}

	if [ {HOST} == 'biowulf.nih.gov' ]; then
		echo "mpiexec  -envall -host `scontrol show hostname ${{SLURM_NODELIST}} | paste -d',' -s` -np ${{SLURM_NTASKS}} novoalignMPI -d {input.index} -f {input.R[0]} {input.R[1]} -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -c 30 -e 100 -F STDFQ --hlimit 7 -i 200 100 -l 30 -o SAM  \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:{params.platform}\" -p 5,2 -t 250  | samtools view -Sbh - | samtools sort -m 30000000000 - {wildcards.subject}/{TIME}/{wildcards.sample}/{wildcards.sample}.novo"

	elif [ {HOST} == 'login01' ]; then
		echo "`cat /cm/local/apps/torque/var/spool/aux/${{PBS_JOBID}} | sort | uniq > ${{LOCAL}}/hosts.txt`
		mpiexec -f ${{LOCAL}}/hosts.txt -np 3 novoalignMPI -d {input.index} -f {input.R[0]} {input.R[1]} -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -c 30 -e 100 -F STDFQ --hlimit 7 -i 200 100 -l 30 -o SAM \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:{params.platform}\"  -p 5,2 -t 250  | samtools view -Sbh - | samtools sort -m 30000000000 - {wildcards.subject}/{TIME}/{wildcards.sample}/{wildcards.sample}.novo"
	else 
		novoalign -c ${{THREADS}} -d {input.index} -f {input.R[0]} {input.R[1]} -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -c 30 -e 100 -F STDFQ --hlimit 7 -i 200 100 -l 30 -o SAM \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:{params.platform}\"  -p 5,2 -t 250  | samtools view -Sbh - | samtools sort -m 30000000000 - {wildcards.subject}/{TIME}/{wildcards.sample}/{wildcards.sample}.novo	
	fi
	novoalign -c ${{THREADS}} -d {input.index} -f {input.R[0]} {input.R[1]} -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -c 30 -e 100 -F STDFQ --hlimit 7 -i 200 100 -l 30 -o SAM \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:{params.platform}\"  -p 5,2 -t 250  | samtools view -Sbh - | samtools sort -m 30000000000 - {wildcards.subject}/{TIME}/{wildcards.sample}/{wildcards.sample}.novo
	samtools index {wildcards.subject}/{TIME}/{wildcards.sample}/{wildcards.sample}.novo.bam
	#######################
	"""
############
#       GenotypeFile
############
# Using Older version of samtools for this purpose
rule GENOTYPING:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.{aligner}.bam",
		interval=config["genotypeBed"],
		ref=config["reference"],
		vcf2genotype=NGS_PIPELINE + "/scripts/vcf2genotype.pl",
		vcf2loh=NGS_PIPELINE + "/scripts/vcf2loh.pl",
	output:
		vcf="{base}/{TIME}/{sample}/calls/{sample}.{aligner}.samtools.vcf",
		gt= "{base}/{TIME}/{sample}/qc/{sample}.{aligner}.gt",
		loh="{base}/{TIME}/{sample}/qc/{sample}.{aligner}.loh"
	version: config["samtools_old"]
	params:
		rulename  = "genotype",
		batch     = config[config['host']]["job_genotype"],
		dest	  = config["genotypeDest"],
		host	  = config["host"]
	shell: """
	#######################
	module load samtools/{version}
	samtools mpileup -u -C50 -f {input.ref} -l {input.interval} {input.bam} | bcftools view -gc - >{output.vcf}
	perl {input.vcf2genotype} {output.vcf} >{output.gt}
	
	perl {input.vcf2loh} {output.vcf} >{output.loh}
	
	#######################
	"""
############
# Genotyping On Sample
############
rule SampleGT:
	input:
		gtFiles=lambda wildcards: SUB_GT[wildcards.subject],
		mail =NGS_PIPELINE + "/scripts/tsv2html.sh",
		score=NGS_PIPELINE + "/scripts/scoreGenotyes.pl"
	output:
		"{subject}/{TIME}/qc/{subject}.genotyping.txt",
	params:
		rulename 	= "SampleGT",
		batch    	= config[config['host']]["job_default"],
		mail     	= config["mail"],
		host	 	= config["host"],
		diagnosis	= lambda wildcards: config['Diagnosis'][SUBJECT_TO_SAMPLE[wildcards.subject][0]]
	shell: """
	#######################
	mkdir -p {wildcards.subject}/{TIME}/qc/GT
	mkdir -p {wildcards.subject}/{TIME}/qc/RATIO/
	cp {input.gtFiles} {wildcards.subject}/{TIME}/qc/GT/
	echo Sample >{wildcards.subject}/{TIME}/qc/RATIO/FirstColumn

	for file in {wildcards.subject}/{TIME}/qc/GT/*
	do
		sample=`basename ${{file}} .gt`
		echo ${{sample}} >>{wildcards.subject}/{TIME}/qc/RATIO/FirstColumn
		echo ${{sample}} >>{wildcards.subject}/{TIME}/qc/RATIO/${{sample}}.ratio
		for file2 in {wildcards.subject}/{TIME}/qc/GT/*
		do
			perl {input.score} ${{file}} ${{file2}} >>{wildcards.subject}/{TIME}/qc/RATIO/${{sample}}.ratio
		done
	done
	paste {wildcards.subject}/{TIME}/qc/RATIO/FirstColumn {wildcards.subject}/{TIME}/qc/RATIO/*.ratio >{wildcards.subject}/{TIME}/qc/{wildcards.subject}.genotyping.txt
	rm -rf {wildcards.subject}/{TIME}/qc/GT/ {wildcards.subject}/{TIME}/qc/RATIO/
	sed -i 's/Sample_//g' {wildcards.subject}/{TIME}/qc/{wildcards.subject}.genotyping.txt
	sed -i 's/.bwa//g' {wildcards.subject}/{TIME}/qc/{wildcards.subject}.genotyping.txt
	sed -i 's/.star//g' {wildcards.subject}/{TIME}/qc/{wildcards.subject}.genotyping.txt
	ssh {params.host} "sh {input.mail} --name {wildcards.subject} --diagnosis '{params.diagnosis}' --head {WORK_DIR}/{wildcards.subject}/{TIME}/qc/{wildcards.subject}.genotyping.txt | /usr/bin/mutt -e \\\"my_hdr Content-Type: text/html\\\" -s 'Genotyping Result on {wildcards.subject}' `whoami`@mail.nih.gov {params.mail} "
	#######################
	"""
############
#       BamQC
############
rule BamQC:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.bwa.final.bam",
		bai="{base}/{TIME}/{sample}/{sample}.bwa.final.bam.bai",
		interval= lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]].replace('.bed', '.gff'),
	output:
		"{base}/{TIME}/{sample}/qc/BamQC/qualimapReport.html"
	version: config["qualimap"]
	params:
		rulename = "bamqc",
		batch	 = config[config['host']]["job_qualimap"],
		outdir	 ="{base}/{TIME}/{sample}/qc/BamQC",
	shell: """
	#######################
	module load qualimap/{version}
	qualimap bamqc -c -bam {input.bam} -outdir {params.outdir} -gff {input.interval} -nt ${{THREADS}} --java-mem-size=${{MEM}}G
	#######################
	"""
############
#       QC_Sum
############
rule QC_Sum:
	input:
		ALL_QC,
		ALL_FASTQC,
		convertor = NGS_PIPELINE + "/scripts/makeQC.pl"
	output:
		"QC_AllSamples.txt"
	version: "v1.1"
	params:
		rulename = "qc_sum",
		batch    = config[config['host']]['job_default']
	shell: """
	#######################
	perl {input.convertor} `pwd` >{output}
	#######################
	"""
############
#       Samtools flagstat
############
rule FLAGSTAT:
	input:	"{base}/{TIME}/{sample}/{sample}.{aligner}.final.bam"
	output: "{base}/{TIME}/{sample}/qc/{sample}.{aligner}.flagstat.txt"
	version: config["samtools"]
	params:
		rulename  = "flagstat",
		batch     = config[config['host']]["job_flagstat"]
	shell: """
	#######################
	module load samtools/{version}
	samtools flagstat {input} > {output}
	#######################
	"""
############
#       Reads Count for every Bed Region on bam file
#	 capture region corrected
############
rule CopyNumber:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.bwa.final.bam",
		interval= lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]].replace("target","targetbp"),
		flagstat="{base}/{TIME}/{sample}/qc/{sample}.bwa.flagstat.txt",
		tool=NGS_PIPELINE+ "/scripts/copyNumber.sh"
	output:
		"{base}/{TIME}/{sample}/copyNumber/{sample}.count.txt"
	version: config["samtools"]
	params:
		rulename  = "CN",
		batch     = config[config['host']]["job_copynumber"]
	shell: """
	#######################
	module load samtools/{version}
	TotalReads=`samtools view -bh -L {input.interval} {input.bam} |samtools flagstat - |head -1|sed -e 's/\s/\\t/g' |cut -f 1`
	#TotalReads=`head -1 {input.flagstat} | sed -e 's/\s/\\t/g' |cut -f 1`
	split -a 5 -d -l 12000 {input.interval} ${{LOCAL}}/input
	for file in ${{LOCAL}}/input*
	do
		sh {input.tool} ${{TotalReads}} ${{file}} {input.bam} ${{file}}.out &
	done
	wait;
	cat ${{LOCAL}}/input*.out >{output}
	#######################
	"""
############
#       Somatic Copy Number
############
rule CN_LRR_old:
	input:
		files=lambda wildcards: somaticCopy[wildcards.Tumor],
		ref=config["gene_coord"],
		index=config["reference"].replace('.fasta', '.index.txt'),
		tool=NGS_PIPELINE+ "/scripts/AddGene.pl",
		cgc    = config["annovar_data"]+"geneLists/combinedList_04292016",
		filter=NGS_PIPELINE+ "/scripts/filterCNV.pl"
	output:
		out=     "{subject}/{TIME}/{Tumor}/copyNumber/{Tumor}.copyNumber.v1.txt",
		hq=      "{subject}/{TIME}/{Tumor}/copyNumber/{Tumor}.hq.v1.txt",
		final=   "{subject}/{TIME}/{Tumor}/copyNumber/{Tumor}.CN.v1.annotated.txt",
		filtered="{subject}/{TIME}/{Tumor}/copyNumber/{Tumor}.CN.v1.filtered.txt"
	params:
		rulename = "LRR",
		batch    = config[config['host']]["job_default"],
		tool     = NGS_PIPELINE+ "/scripts/ListStatistics.R"
	shell: """
	#######################
	module load R
	module load bedtools/2.25.0
	mkdir -p {wildcards.subject}/{TIME}/Actionable/
	echo -e "#Chr\\tStart\\tEnd\\tNormalCoverage\\tTumorCoverage\\tRatio\\tLRR\\tGene(s)\\tStrand(s)" >{output.out}
	paste {input.files} |cut -f 1-4,8 |awk '{{OFS="\\t"}};{{print $1,$2,$3,$4,$5,($5+1)/($4+1),log(($5+1)/($4+1))/log(2)}}' >{output.out}.temp1

	intersectBed -a {input.ref} -b {input.files[0]} >{output.out}.temp
	perl {input.tool} {output.out}.temp {output.out}.temp1 >>{output.out}

	awk '{{if($4>=30) print $0}}' {output.out} >{output.hq}
	median=`cut -f7 {output.hq}|grep -v LRR | {params.tool} --stat median --file - `
	MAD=`cut -f7 {output.hq}|grep -v LRR | {params.tool} --stat MAD --file - `
	min=`echo "${{median}}-(2*${{MAD}})"|bc`
	max=`echo "${{median}}+(2*${{MAD}})"|bc`

	perl {input.filter} filter {output.out} ${{min}} ${{max}} {input.cgc} |sortBed -faidx {input.index} -header -i - >{output.final}
	cp -f {output.final} {wildcards.subject}/{TIME}{ACT_DIR}{wildcards.Tumor}.copyNumber.v1.txt
	geneList=`grep -P "Gain|Loss" {output.final} |cut -f 11 |sort |uniq |grep -v "^-$"`

	head -1 {output.final} >{output.filtered}
	for gene in ${{geneList}};
	do
		awk -v gene=${{gene}} '{{if($11 == gene) print $0}}' {output.final}
	done |sort |uniq |sortBed -faidx {input.index} -header -i - >>{output.filtered}
	cp -f {output.filtered} {wildcards.subject}/{TIME}{ACT_DIR}{wildcards.Tumor}.CN.v1.filtered.txt
	rm -rf {output.out}.temp1 {output.out}.temp
	#######################
	"""
############
#       Somatic Copy Number LRR (Median Corrected)
############
rule CN_LRR:
	input:
		files=lambda wildcards: somaticCopy[wildcards.Tumor],
		ref=config["gene_coord"],
		index=config["reference"].replace('.fasta', '.index.txt'),
		tool=NGS_PIPELINE+ "/scripts/AddGene.pl",
		cgc    = config["annovar_data"]+"geneLists/combinedList_04292016",
		filter=NGS_PIPELINE+ "/scripts/filterCNV.pl"
	output:
		out=     "{subject}/{TIME}/{Tumor}/copyNumber/{Tumor}.copyNumber.txt",
		hq=      "{subject}/{TIME}/{Tumor}/copyNumber/{Tumor}.hq.txt",
		final=   "{subject}/{TIME}/{Tumor}/copyNumber/{Tumor}.CN.annotated.txt",
		filtered="{subject}/{TIME}/{Tumor}/copyNumber/{Tumor}.CN.filtered.txt"
	params:
		rulename = "LRR",
		batch    = config[config['host']]["job_default"],
		tool     = NGS_PIPELINE+ "/scripts/ListStatistics.R"
	shell: """
	#######################
	module load R
	module load bedtools/2.25.0
	mkdir -p {wildcards.subject}/{TIME}/Actionable/
	echo -e "#Chr\\tStart\\tEnd\\tNormalCoverage\\tTumorCoverage\\tRatio\\tLRR\\tGene(s)\\tStrand(s)" >{output.out}
	paste {input.files} |cut -f 1-4,8 |awk '{{OFS="\\t"}};{{print $1,$2,$3,$4,$5,($5+1)/($4+1),log(($5+1)/($4+1))/log(2)}}' >{output.out}.temp1

	intersectBed -a {input.ref} -b {input.files[0]} >{output.out}.temp
	perl {input.tool} {output.out}.temp {output.out}.temp1 >>{output.out}

	# Created the first file

	awk '{{if($4>=30) print $0}}' {output.out} >{output.hq}
	median=`cut -f7 {output.hq}|grep -v LRR | {params.tool} --stat median --file - `

	corr_factor=`echo "0 - ${{median}}"|bc`
	echo -e "#Chr\\tStart\\tEnd\\tNormalCoverage\\tTumorCoverage\\tRatio\\tLRR\\tGene(s)\\tStrand(s)" >{output.out}.corrected
	grep -v Ratio {output.out} |awk -v factor=${{corr_factor}} '{{OFS="\\t"}};{{print $1,$2,$3,$4,$5,$6,$7 + (factor),$8,$9}}' >>{output.out}.corrected

	# Created corrected file
	awk '{{if($4>=30) print $0}}' {output.out}.corrected >{output.hq}

	median=`cut -f7 {output.hq}| grep -v LRR |sort -n |{params.tool} --stat median --file - `

	MAD=`cut -f7 {output.hq}|grep -v LRR |{params.tool} --stat MAD --file - `
	min=`echo "${{median}}-(2.5*${{MAD}})"|bc`
	max=`echo "${{median}}+(2.5*${{MAD}})"|bc`

	perl {input.filter} filter {output.out}.corrected ${{min}} ${{max}} {input.cgc} |sortBed -faidx {input.index} -header -i - >{output.final}
	cp -f {output.final} {wildcards.subject}/{TIME}{ACT_DIR}{wildcards.Tumor}.copyNumber.txt
	set +eo pipefail
	geneList=`grep -P "Gain|Loss" {output.final} |cut -f 11 |sort |uniq |grep -v "^-$"`

	head -1 {output.final} >{output.filtered}
	if [ -n "${{geneList}}" ]; then
		for gene in ${{geneList}};
		do
			awk -v gene=${{gene}} '{{if($11 == gene) print $0}}' {output.final}
		done |sort |uniq |sortBed -faidx {input.index} -header -i - >>{output.filtered}
		cp -f {output.filtered} {wildcards.subject}/{TIME}{ACT_DIR}{wildcards.Tumor}.CN.v2.filtered.txt
		rm -rf {output.out}.temp1 {output.out}.temp
	else
		echo "No Gains or Losses Found in this sample" >>{output.filtered}
		cp -f {output.filtered} {wildcards.subject}/{TIME}{ACT_DIR}{wildcards.Tumor}.CN.v2.filtered.txt
                rm -rf {output.out}.temp1 {output.out}.temp
	fi
	#######################
	"""
############
#       Hotspot Coverage
############
rule HotSpotCoverage:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.{aligner}.final.bam",
		interval=config["hotspot_intervals"]
	output: "{base}/{TIME}/{sample}/qc/{sample}.{aligner}.hotspot.depth"
	version: config["bedtools"]
	params:
		rulename  = "HotSpotCov",
		batch     = config[config['host']]["job_hotspot"],
		samtools  = config["samtools"]
	shell: """
	#######################
	module load samtools/{params.samtools}
	module load bedtools/{version}
	samtools view -hF 0x400 -q 30 {input.bam} | samtools view -ShF 0x4 - | samtools view -SuF 0x200 - | bedtools coverage -abam - -b {input.interval} >{output}
	#######################
	"""
############
# Coverage
############
rule Coverage:
	input:
		bam="{subject}/{TIME}/{sample}/{sample}.bwa.final.bam",
		bai="{subject}/{TIME}/{sample}/{sample}.bwa.final.bam.bai",
		interval= lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]],
	output:
		"{subject}/{TIME}/{sample}/qc/{sample}.bwa.coverage.txt"
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
# CoveragePlot
############
rule CoveragePlot:
	input:
		covFiles=lambda wildcards: SUB_COV[wildcards.subject],
		coverage =NGS_PIPELINE + "/scripts/coverage.R"
	output: "{subject}/{TIME}/qc/{subject}.coveragePlot.png",
	version: config["R"]
	params:
		rulename = "covplot",
		batch    = config[config['host']]["job_covplot"]
	shell: """
	#######################

	cp -f {input.covFiles} ${{LOCAL}}

	module load R
	R --vanilla --slave --silent --args ${{LOCAL}} {output} {wildcards.subject} <{input.coverage}
	#######################
	"""
############
# IGV Session file
############
rule IGV_Session:
	input: bams=lambda wildcards: SUB_IGV[wildcards.subject]
	output: "{subject}/{TIME}/igv/session_{subject}.xml"
	message: "Making IGV session xml file for {wildcards.subject}"
	params:
		rulename = "igv_session",
		batch    = config[config['host']]["job_covplot"],
		work_dir =  WORK_DIR
	shell: """
	#######################
	dir=`echo {params.work_dir} | sed -e 's/\/data\/khanlab/K:/g' | sed -e 's/\/projects\/Clinomics/Y:/g'`
	echo "<?xml version=\\"1.0\\" encoding=\\"UTF-8\\"?>" >{output}
	echo "<Global genome=\\"hg19\\" locus=\\"\\" version=\\"3\\">" >>{output}
	echo "\t<Resources>" >>{output}
	for BAM in {input.bams}
	do
		bam=`echo "${{dir}}/${{BAM}}" |sed -e 's/\//\\\\\\/g'`
		echo "\t\t<Resource path=\\"${{bam}}\\"/>" >>{output}
	done
	echo "\t</Resources>" >>{output}
	echo "</Global>" >>{output}
	#######################
	"""
############
# Circos Plot
############
rule Circos:
	input:
		lohFiles=lambda wildcards: SUB_LOH[wildcards.subject],
		circos =NGS_PIPELINE + "/scripts/circos.R"
	output:
		"{subject}/{TIME}/qc/{subject}.circos.png",
	version: config["R"]
	params:
		rulename = "Circos",
		batch    = config[config['host']]["job_covplot"]
	shell: """
	#######################
	cp -f {input.lohFiles} ${{LOCAL}}
	module load R
	R --vanilla --slave --silent --args ${{LOCAL}} {output} {wildcards.subject} <{input.circos}
	#######################
	"""
############
# Box Plot Hotspot
############
rule BoxPlot_Hotspot:
	input:
		covFiles=lambda wildcards: SUB_HOT[wildcards.subject],
		boxplot =NGS_PIPELINE + "/scripts/boxplot.R"
	output:
		"{subject}/{TIME}/qc/{subject}.hotspot_coverage.png",
	version: config["R"]
	params:
		rulename = "Boxplot",
		batch    = config[config['host']]["job_covplot"]
	shell: """
	#######################
	cp -f {input.covFiles} ${{LOCAL}}
	module load R
	R --vanilla --slave --silent --args ${{LOCAL}} {output} {wildcards.subject} <{input.boxplot}
	#######################
	"""
############
#       Picard Mark Duplicates
############
rule Picard_MarkDup:
	input:
		bam="{subject}/{TIME}/{sample}/{sample}.{base}.bam",
		bai="{subject}/{TIME}/{sample}/{sample}.{base}.bam.bai"

	output:
		bam=temp("{subject}/{TIME}/{sample}/{sample}.{base}.dd.bam"),
		index=temp("{subject}/{TIME}/{sample}/{sample}.{base}.dd.bam.bai"),
		metrics="{subject}/{TIME}/{sample}/qc/{base}.markdup.txt"
	version: config["picard"]
	params:
		rulename  = "mark_dup",
		batch     = config[config['host']]["job_markdup"],
		samtools  = config["samtools"]
	shell: """
	#######################
	module load picard/{version}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $PICARD_JAR MarkDuplicates AS=true M={output.metrics} I={input.bam} O={output.bam} REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT
	module load samtools/{params.samtools}
	samtools index {output.bam}
	######################
	"""
############
# copy novo bam to novo.final bam
############
rule CopyNovoBam:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.novo.dd.bam",
		bai="{base}/{TIME}/{sample}/{sample}.novo.dd.bam.bai",
	output:
		bam="{base}/{TIME}/{sample}/{sample}.novo.final.bam",
		bai="{base}/{TIME}/{sample}/{sample}.novo.final.bam.bai",
	params:
		rulename = "copy",
		batch    = config[config['host']]['job_default']
	shell: """
	#######################

	mv {input.bam} {output.bam}
	mv {input.bai} {output.bai}

	#######################
	"""
############
#       GATK Best Practices
############
rule GATK:
	input: 	bam="{base}/{TIME}/{sample}/{sample}.bwa.dd.bam",
		bai="{base}/{TIME}/{sample}/{sample}.bwa.dd.bam.bai",
		ref=config["reference"],
		phase1=config["1000G_phase1"],
		mills=config["Mills_and_1000G"]
	output:
		bam="{base}/{TIME}/{sample}/{sample}.bwa.final.bam",
		index="{base}/{TIME}/{sample}/{sample}.bwa.final.bam.bai",
	version: config["GATK"]
	params:
		rulename  = "gatk",
		batch     = config[config['host']]["job_gatk"]
	shell: """
	#######################
	module load GATK/{version}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T RealignerTargetCreator -R {input.ref} -known {input.phase1} -known {input.mills} -I {input.bam} -o ${{LOCAL}}/{wildcards.sample}.realignment.intervals
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T IndelRealigner -R {input.ref} -known {input.phase1} -known {input.mills} -I {input.bam} --targetIntervals ${{LOCAL}}/{wildcards.sample}.realignment.intervals -o ${{LOCAL}}/{wildcards.sample}.lr.bam
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T BaseRecalibrator -R {input.ref} -knownSites {input.phase1} -knownSites {input.mills} -I ${{LOCAL}}/{wildcards.sample}.lr.bam -o ${{LOCAL}}/{wildcards.sample}.recalibration.matrix.txt
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T PrintReads -R {input.ref} -I ${{LOCAL}}/{wildcards.sample}.lr.bam -o {output.bam} -BQSR ${{LOCAL}}/{wildcards.sample}.recalibration.matrix.txt
	mv -f {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}.bwa.final.bai {output.index}
	######################
	"""
############
#	Bam2MPG
############
rule Bam2MPG:
	input:
		bam="{subject}/{TIME}/{sample}/{sample}.novo.final.bam",
		bai="{subject}/{TIME}/{sample}/{sample}.novo.final.bam.bai",
		ref=config["reference"],
		interval=lambda wildcards: config['target_intervals'][config['sample_captures'][wildcards.sample]].replace("target","targetbp")
	output:
		snps="{subject}/{TIME}/{sample}/calls/{sample}.bam2mpg.vcf.gz",
		vcf="{subject}/{TIME}/{sample}/calls/{sample}.bam2mpg.raw.vcf"
	version: config["bam2mpg"]
	params:
		rulename  = "bam2mpg",
		batch     = config[config['host']]["job_bam2mpg"],
		samtools  = config["samtools"],
		vcftools  = config["vcftools"],
		parallel  = NGS_PIPELINE + "/scripts/parallel"
	shell: """
	#######################
	if [ -f {wildcards.subject}/{TIME}/{wildcards.sample}/calls/{wildcards.sample}.bam2mpg.vcf.gz ]; then
        	rm -rf {wildcards.subject}/{TIME}/{wildcards.sample}/calls/{wildcards.sample}.bam2mpg.vcf.*

	fi

	module load bam2mpg/{version}
	module load vcftools/{params.vcftools}
	for CHR in `seq 1 22` X Y;
	do
	echo "bam2mpg --qual_filter 20 -bam_filter '-q31' --region chr${{CHR}} --only_nonref --snv_vcf ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.snps.vcf --div_vcf ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.indel.vcf {input.ref} {input.bam}"
	done >{wildcards.subject}/{TIME}/{wildcards.sample}/calls/swarm.file
	cat   {wildcards.subject}/{TIME}/{wildcards.sample}/calls/swarm.file |{params.parallel} -j 10 --no-notice
	rm -rf {wildcards.subject}/{TIME}/{wildcards.sample}/calls/swarm.file
	
	for CHR in `seq 1 22` X Y
	do
		
		echo "Started indexing chr${{CHR}}"
		cat ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.snps.vcf | vcf-sort >${{LOCAL}}/chr${{CHR}}{wildcards.sample}.snps.tmp.vcf
		mv -f ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.snps.tmp.vcf ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.snps.vcf
		bgzip ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.snps.vcf
		cat ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.indel.vcf | vcf-sort >${{LOCAL}}/chr${{CHR}}{wildcards.sample}.indel.tmp.vcf
		mv -f ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.indel.tmp.vcf ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.indel.vcf
		bgzip ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.indel.vcf
		tabix -p vcf ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.snps.vcf.gz
		tabix -p vcf ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.indel.vcf.gz
		echo "Finished indexing chr${{CHR}}"
	done
	echo "Combine chr level vcf files"
	vcf-concat ${{LOCAL}}/chr*{wildcards.sample}.*.vcf.gz >${{LOCAL}}/{wildcards.sample}.snps.vcf
	echo "Restrict to Bed file"

	vcftools --vcf ${{LOCAL}}/{wildcards.sample}.snps.vcf --bed {input.interval} --out ${{LOCAL}}/{wildcards.sample} --recode --keep-INFO-all

	sed -e 's/SAMPLE/{wildcards.sample}/g' ${{LOCAL}}/{wildcards.sample}.recode.vcf |vcf-sort >{wildcards.subject}/{TIME}/{wildcards.sample}/calls/{wildcards.sample}.bam2mpg.vcf

	bgzip {wildcards.subject}/{TIME}/{wildcards.sample}/calls/{wildcards.sample}.bam2mpg.vcf
	tabix -f -p vcf {wildcards.subject}/{TIME}/{wildcards.sample}/calls/{wildcards.sample}.bam2mpg.vcf.gz
	gunzip -c {wildcards.subject}/{TIME}/{wildcards.sample}/calls/{wildcards.sample}.bam2mpg.vcf.gz >{wildcards.subject}/{TIME}/{wildcards.sample}/calls/{wildcards.sample}.bam2mpg.raw.vcf
	#######################
	"""
############
# Subject Bam2MPG
############
rule Sub_MPG:
	input:
		vcf=lambda wildcards: SUB_MPG[wildcards.subject],
	output:
		vcf="{subject}/{TIME}/{subject}/calls/{subject}.bam2mpg.raw.vcf"
	version: config["vcftools"]
	params:
		rulename = "mergevcf",
		batch    = config[config['host']]["job_default"],
		vcforder = NGS_PIPELINE + "/scripts/vcfOrderCol.R"
	shell: """
	#######################
	module load vcftools/{version}
	module load R
	vcf-merge {input.vcf} > {output.vcf}.tmp
	{params.vcforder} -i {output.vcf}.tmp -o {output.vcf}
	rm -rf {output.vcf}.tmp
	#######################
	"""
############
#       MuTect
############
rule MuTect:
	input:
		lambda wildcards: somaticPairs[wildcards.Tumor],
		ref=config["reference"],
		dbsnp=config["dbsnp"],
		cosmic=config["cosmic"],
		interval=lambda wildcards: config['target_intervals'][pairedCapture[wildcards.Tumor]].replace("target","targetbp")
	output:
		vcf="{subject}/{TIME}/{Tumor}/calls/{Tumor}.MuTect.raw.vcf",
		call_stats="{subject}/{TIME}/{Tumor}/qc/{Tumor}.mutect.call_stats.txt",
		coverage="{subject}/{TIME}/{Tumor}/qc/{Tumor}.mutect.coverage.wig.txt"
	version: config["MuTect"]
	params:
		rulename = "MuTect",
		batch    = config[config['host']]["job_mutect"],
		vcforder = NGS_PIPELINE + "/scripts/vcfOrderCol.R",
		mt       = "--max_alt_allele_in_normal_fraction 0.05 --max_alt_alleles_in_normal_count 4 --min_qscore 20 -rf MappingQuality -mmq 30"
	shell: """
	#######################
	module load muTect/{version}
	module load R
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $MUTECT_JAR -T MuTect \
		--reference_sequence {input.ref} \
		--cosmic {input.cosmic} \
		--dbsnp {input.dbsnp} \
		--input_file:normal {input[2]} \
		--input_file:tumor {input[0]} \
		--intervals  {input.interval} \
		--coverage_file {output.coverage} \
		--out {output.call_stats} \
		--vcf {output.vcf}.vcf \
		{params.mt}

	{params.vcforder} -i {output.vcf}.vcf -o {output.vcf}
	rm -rf {output.vcf}.vcf
	#######################
	"""
############
#       Strelka
############
rule Strelka:
	input:
		lambda wildcards: somaticPairs[wildcards.Tumor],
		ref=config["reference"],
		config=config["strelka_config"],
		interval=lambda wildcards: config['target_intervals'][pairedCapture[wildcards.Tumor]].replace("target","targetbp")
	output:
		snps="{subject}/{TIME}/{Tumor}/calls/{Tumor}.strelka.snvs.raw.vcf",
		indels="{subject}/{TIME}/{Tumor}/calls/{Tumor}.strelka.indels.raw.vcf"
	version: config["strelka"]
	params:
		rulename = "Strelka",
		batch    = config[config['host']]["job_strelka"],
		vcftools = config["vcftools"]
	shell: """
	#######################
	module load strelka/{version}
	configureStrelkaWorkflow.pl --normal={input[2]} --tumor={input[0]}\
	--ref={input.ref} --config={input.config} --output-dir=${{LOCAL}}/strelka
	make -j ${{SLURM_CPUS_PER_TASK}} -f ${{LOCAL}}/strelka/Makefile
	module load vcftools/{params.vcftools}
	vcftools --vcf ${{LOCAL}}/strelka/results/passed.somatic.snvs.vcf --bed {input.interval} --out {output.snps} --recode --keep-INFO-all
	mv -f {output.snps}.recode.vcf {output.snps}
	vcftools --vcf ${{LOCAL}}/strelka/results/passed.somatic.indels.vcf --bed {input.interval}  --out {output.indels} --recode --keep-INFO-all
	mv -f {output.indels}.recode.vcf {output.indels}
	NORMAL=`basename {input[2]} .bwa.final.bam`
	sed -i "s/FORMAT\\tNORMAL\\tTUMOR/FORMAT\\t${{NORMAL}}\\t{wildcards.Tumor}/g" {output.snps}
	sed -i "s/FORMAT\\tNORMAL\\tTUMOR/FORMAT\\t${{NORMAL}}\\t{wildcards.Tumor}/g" {output.indels}

	#######################
	"""
############
# Subject Hapcaller
############
rule Sub_HapCall:
	input:
		bams=lambda wildcards: SUB_BAMS[wildcards.subject],
		ref=config["reference"],
		dbsnp=config["dbsnp"],
		interval= config["coding_bed"]
	output:
		vcf="{subject}/{TIME}/{subject}/calls/{subject}.hapCaller.raw.vcf"
	version: config["GATK"]
	params:
		rulename = "HC",
		batch    = config[config['host']]["job_gatk"]
	run:
		bamFiles = " ".join('-I ' + bam for bam in input.bams)
		shell("""
	#######################
	module load GATK/{version}
	gawk '{{print $1 "\t" $2-1 "\t" $3}}' {input.interval} > ${{LOCAL}}/target_intervals.bed
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T HaplotypeCaller -R {input.ref} {bamFiles} -L ${{LOCAL}}/target_intervals.bed -o {output.vcf} --dbsnp {input.dbsnp} -mbq 20 -mmq 30
	#######################
	""")
############
# Subject Platypus
############
rule Sub_Platypus:
	input:
		bams=lambda wildcards: SUB_BAMS[wildcards.subject],
		ref=config["reference"],
		dbsnp=config["dbsnp"],
		interval=config["coding_bed"]
	output:
		vcf="{subject}/{TIME}/{subject}/calls/{subject}.platypus.raw.vcf"
	version: config["platypus"]
	log: "log/platypus.{subject}"
	params:
		rulename = "PLAT",
		batch    = config[config['host']]["job_platypus"]
	shell: """
	#######################
	module load platypus/{version}
	LIST=`echo {input.bams}|sed -e 's/ /,/g'`
	platypus callVariants --nCPU=${{THREADS}} --bufferSize=1000000 --maxReads=100000000 --bamFiles=${{LIST}} --regions={input.interval} --output={output.vcf} --refFile={input.ref}  --logFileName={log}
	sed -i 's/.bwa.final//g' {output.vcf}
	#######################
	"""
############
#	snpEff
############
rule SNPEff:
	input:
		vcf="{subject}/{TIME}/calls/{base}.raw.vcf",
		ref=config["reference"],
		snpEff_config=config["snpEff_config"],
	output:
		eff="{subject}/{TIME}/calls/{base}.raw.snpEff.vcf"
	version: config["snpEff"]
	params:
		rulename      ="snpEff",
		batch	      =config[config['host']]["job_snpeff"],
		snpEff_genome =config["snpEff_genome"],
		annovar       =config["annovar"]
	shell: """
	#######################
	module load snpEff/{version}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $SNPEFF_JARPATH/SnpSift.jar dbnsfp -c {input.snpEff_config} -a {input.vcf} | java -Xmx${{MEM}}g -jar $SNPEFF_JARPATH/snpEff.jar -t -canon {params.snpEff_genome} > {output.eff}
	#######################
	"""
############
#       vcf2txt
############
rule VCF2TXT:
	input:
		eff="{subject}/{TIME}/calls/{base}.raw.snpEff.vcf",
		vcf2txt=NGS_PIPELINE + "/scripts/vcf2txt.pl"
	output:
		txt="{subject}/{TIME}/calls/{base}.snpEff.txt"
	params:
		rulename      ="vcf2txt",
		batch         =config[config['host']]["job_default"],
		annovar       =config["annovar"]
	shell: """
	#######################
	module load annovar/{params.annovar}
	perl {input.vcf2txt} {input.eff} ${{LOCAL}} >{output.txt}
	#######################
	"""

############
#	MakeList
############
rule FormatInput:
	input:
		txtFiles=lambda wildcards: SUBJECT_VCFS[wildcards.subject],
		convertor= NGS_PIPELINE + "/scripts/MakeAnnotationInputs.pl"
	output:
		temp("{subject}/annotation/AnnotationInput.anno"),
		temp("{subject}/annotation/AnnotationInput.sift")
	version: config["annovar"]
	params:
		rulename   = "FormatInput",
		batch      = config[config['host']]["job_default"],
		fAEV       = NGS_PIPELINE + "/scripts/findAlreadyExistingVariants.pl"
	shell: """
	#######################
	module load annovar/{version}
	cut -f 1-5 {input.txtFiles} |sort |uniq > {wildcards.subject}/annotation/AnnotationInput
	#cut -f 1-5 {input.txtFiles} |sort |uniq > {wildcards.subject}/annotation/allSites
	#perl {params.fAEV} {wildcards.subject}/annotation/allSites annovar/AnnotationInput.final.txt {wildcards.subject}/annotation/AnnotationInput.annotations.final.txt {wildcards.subject}/annotation/AnnotationInput
	perl {input.convertor} {wildcards.subject}/annotation/AnnotationInput
	rm -rf "{wildcards.subject}/annotation/AnnotationInput.pph",
	#######################
	"""
############
#	Custom Annotation
############
rule Annotation:
	input:
		"{subject}/annotation/AnnotationInput.anno",
		TableAnnovar=NGS_PIPELINE + "/scripts/TableAnno.sh",
		custom     =NGS_PIPELINE + "/scripts/addAnnotation.pl"
	output:
		temp("{subject}/annotation/AnnotationInput.docm")
	version: config["annovar"]
	params:
		rulename   = "Annotation",
		batch      = config[config['host']]["job_annovar"],
		RefData    = config["annovar_data"],
		build      = config["build"],
	shell: """
	#######################
	module load annovar/{version}
	sh {input.TableAnnovar} {wildcards.subject}/annotation AnnotationInput {input.custom} {params.RefData}
	#######################
	"""
############
#       SIFT
############
rule SIFT:
	input:
		sift="{subject}/annotation/AnnotationInput.sift",
		convertor  = NGS_PIPELINE + "/scripts/ParseSIFT.pl"
	output:
		temp("{subject}/annotation/AnnotationInput.sift.out")
	version: config["SIFT"]
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
		SIFT_exome_nssnvs.pl -i ${{FILE}} -d $SIFTDB/Human_db_37 -o ${{LOCAL}}  -z ${{DIR}}/`dirname {input.sift}`/${{FILE}}.sift_predictions.tsv
		perl {input.convertor} ${{DIR}}/`dirname {input.sift}`/${{FILE}}.sift_predictions.tsv >${{DIR}}/`dirname {input.sift}`/${{FILE}}.out
	else
		echo -e "Chr\\tStart\\tEnd\\tRef\\tAlt\\tSIFT Prediction\\tSIFT Score" >{output}
	fi
	#######################
	"""
############
#       PPH2   ** Not in use**
############
rule PPH2:
	input:
		pph="{subject}/annotation/AnnotationInput.pph",
		convertor  = NGS_PIPELINE + "/scripts/ParsePPH2.pl"
	output: "{subject}/annotation/AnnotationInput.pph2.out"
	version: config["polyphen2"]
	params:
		rulename   = "PPH2",
		batch      = config[config['host']]["job_PPH2"],
	shell: """
	#######################
	module load polyphen2/{version}
	if [ -s {input.pph} ]; then
		mapsnps.pl -c -g hg19 -U -y {input.pph}.intermediate {input.pph}
		pph_swarm.pl {input.pph}.intermediate -d /scratch/`whoami`/${{RANDOM}}${{RANDOM}} -o {wildcards.subject}/annotation/AnnotationInput.pph2.intermediate.txt --partition ${{SLURM_JOB_PARTITION}} --block
		perl {input.convertor}  {wildcards.subject}/annotation/AnnotationInput.pph2.intermediate.txt >{output}
	else
		touch {output}
	fi
	rm -rf {wildcards.subject}/annotation/AnnotationInput.pph.inter*
	#######################
	"""
############
#	Combine Annotation
############
rule CombineAnnotation:
	input:
		anno="{subject}/annotation/AnnotationInput.docm",
		sift="{subject}/annotation/AnnotationInput.sift.out",
		convertor  = NGS_PIPELINE + "/scripts/CombineAnnotations.pl",
		geneanno   = NGS_PIPELINE + "/scripts/GeneAnnotation.pl",
		filter     = NGS_PIPELINE + "/scripts/filterVariants.pl",
		coding     = NGS_PIPELINE + "/scripts/ProteinCoding.pl",
		blacklisted      = config["annovar_data"]+ "hg19_blacklistedSites.txt"
	output: "{subject}/annotation/{subject}.Annotations.coding.rare.txt"
	version: "1.0"
	params:
		rulename   = "combine",
		batch	   = config[config['host']]["job_Combine"],
		dataDir    = config["annovar_data"]
	shell: """
	#######################
	echo "{wildcards.subject}/annotation/AnnotationInput
{wildcards.subject}/annotation/AnnotationInput.anno.gene
{wildcards.subject}/annotation/AnnotationInput.anno.exac.3
{wildcards.subject}/annotation/AnnotationInput.anno.clinseq
{wildcards.subject}/annotation/AnnotationInput.anno.cadd
{wildcards.subject}/annotation/AnnotationInput.sift.out
{wildcards.subject}/annotation/AnnotationInput.clinvar
{wildcards.subject}/annotation/AnnotationInput.anno.cosmic
{wildcards.subject}/annotation/AnnotationInput.hgmd
{wildcards.subject}/annotation/AnnotationInput.match
{wildcards.subject}/annotation/AnnotationInput.docm
{wildcards.subject}/annotation/AnnotationInput.candl
{wildcards.subject}/annotation/AnnotationInput.tcc
{wildcards.subject}/annotation/AnnotationInput.mcg
{wildcards.subject}/annotation/AnnotationInput.civic
{wildcards.subject}/annotation/AnnotationInput.anno.pcg" >{wildcards.subject}/annotation/list
	perl {input.convertor} {wildcards.subject}/annotation/list >{output}
	perl {input.geneanno} {params.dataDir}hg19_ACMG.txt {output} >>{wildcards.subject}/annotation/AnnotationInput.annotations.final.txt
	perl {input.coding} {wildcards.subject}/annotation/AnnotationInput.annotations.final.txt | perl {input.filter} - {input.blacklisted} 0.05 |sort -n |uniq >{output}.all
	grep -P "Chr\\tStart\\tEnd\\tRef\\tAlt" {output}.all >{output}
	grep -v -P "Chr\\tStart\\tEnd\\tRef\\tAlt" {output}.all >>{output}
	rm -rf {output}.all {wildcards.subject}/annotation/list
	

	rm -rf {wildcards.subject}/annotation/AnnotationInput.pph {wildcards.subject}/annotation/AnnotationInput.anno.* {wildcards.subject}/annotation/AnnotationInput.hgmd {wildcards.subject}/annotation/AnnotationInput.match {wildcards.subject}/annotation/AnnotationInput.candl {wildcards.subject}/annotation/AnnotationInput.tcc {wildcards.subject}/annotation/AnnotationInput.mcg {wildcards.subject}/annotation/AnnotationInput.civic {wildcards.subject}/annotation/AnnotationInput.anno.pcg {wildcards.subject}/annotation/AnnotationInput.clinvar {wildcards.subject}/annotation/AnnotationInput {wildcards.subject}/annotation/allSites
	#######################
	"""
############
#	Add Annotation back to sample level file
############
rule AttachAnnotation:
	input:
		txt="{subject}/{TIME}/{base1}/calls/{base}.snpEff.txt",
		ref="{subject}/annotation/{subject}.Annotations.coding.rare.txt",
		convertor  = NGS_PIPELINE + "/scripts/addAnnotations2vcf.pl"
	output:
		txt="{subject}/{TIME}/{base1}/calls/{base}.annotated.txt"
	version: "1.0"
	params:
		rulename   = "add",
		batch      = config[config['host']]["job_addbackann"],
	shell: """
	#######################
	perl {input.convertor} {input.ref}  {input.txt} >{output.txt}
	#######################
	"""
############
#       Expressed
############
rule Expressed:
	input:
		RNASeq = lambda wildcards: expressedPairs[wildcards.sample],
		Mutation="{subject}/{TIME,[0-9]+}/{sample}/calls/{base}.annotated.txt",
		convertor = NGS_PIPELINE + "/scripts/ExpressedMutations.pl",
		convertor1 = NGS_PIPELINE + "/scripts/mpileup.pl"
	output: "{subject}/{TIME,[0-9]+}/{sample}/calls/{base}.annotated.expressed.txt"
	version: config["samtools"]
	params:
		rulename  = "Expressed",
		batch     = config[config['host']]["job_expressed"],
		name      = lambda wildcards: config["sample_RNASeq"][wildcards.sample]
	shell: """
	#######################
	if [ {output} == 'indels' ]
	then
		perl {input.convertor} {input.Mutation} {input.RNASeq} {params.name} >{output}
	else
		perl {input.convertor1} {input.Mutation} {wildcards.subject}/{TIME}/{params.name}/{params.name}.star.final.bam >{output}	
	fi
	#######################
	"""
############
#       Database Input
############
rule DBinput:
	input:
		txtFiles=lambda wildcards: SUBJECT_ANNO[wildcards.subject][wildcards.group],
		convertor=NGS_PIPELINE + "/scripts/makeDBVariantFile.pl",
		tool=NGS_PIPELINE + "/scripts/AddSampleType.pl"
	output: "{subject}/{TIME,[0-9]+}/{subject}/db/{subject}.{group}"
	params:
		rulename = "makeDBinput",
		batch    = config[config['host']]['job_default'],
		hash 	 = config["sample_type"].items(),
		hash1	 = config["sample_captures"].items(),
	shell: """
	#######################
	perl {input.convertor} {input.txtFiles} |perl {input.tool} - "{params.hash}" "{params.hash1}" >{output}
	#######################
	"""
############
#       Actionable
############
rule Actionable_Somatic:
	input:
		somatic       = "{subject}/{TIME,[0-9]+}/{subject}/db/{subject}.somatic",
		annotation    = "{subject}/annotation/{subject}.Annotations.coding.rare.txt",
		refFile       = config["annovar_data"]+"hg19_SomaticActionableSites.txt",
		cgc           = config["annovar_data"]+"geneLists/CancerGeneCensus.v76.txt",
		combinedList  = config["annovar_data"]+"geneLists/combinedList_04292016",
		annotate      = NGS_PIPELINE + "/scripts/addAnnotations2vcf.pl",
		convertor     = NGS_PIPELINE + "/scripts/" + config["Actionable_mutation"],
	output:
		somatic="{subject}/{TIME,[0-9]+}{ACT_DIR}{subject}.somatic.actionable.txt",
	params:
		rulename  = "ActionableMutations",
		batch    = config[config['host']]['job_default']
	shell: """
	#######################
	perl {input.convertor} somatic  {input.refFile} {input.cgc} {input.combinedList} {input.somatic} {input.annotation} >{output.somatic}
	#######################
	"""
############
#       Actionable
############
rule Actionable_Germline:
	input:
		germline  ="{subject}/{TIME,[0-9]+}/{subject}/db/{subject}.germline",
		annotation="{subject}/annotation/{subject}.Annotations.coding.rare.txt",
		annotate =NGS_PIPELINE + "/scripts/addAnnotations2vcf.pl",
		convertor=NGS_PIPELINE + "/scripts/" + config["Actionable_mutation"],
		cancerGeneCensus = config["annovar_data"]+"geneLists/CGCensus_Hereditary.txt",
		hotspot= config["annovar_data"]+"hg19_SomaticActionableSites.txt",
		combinedList  = config["annovar_data"]+"geneLists/combinedList_04292016",
		combine=NGS_PIPELINE + "/scripts/germlineOnly.pl",
		cgc           = config["annovar_data"]+"geneLists/CancerGeneCensus.v76.txt"
	output:
		germline="{subject}/{TIME,[0-9]+}{ACT_DIR}{subject}.germline.actionable.txt",
	params:
		rulename  = "ActionableMutations",
		batch    = config[config['host']]['job_default']
	shell: """
	#######################
	if [ -e {wildcards.subject}/{TIME}/{wildcards.subject}/db/{wildcards.subject}.somatic ]
	then
		perl {input.convertor} germline {wildcards.subject}/{TIME}/{wildcards.subject}/db/{wildcards.subject}.somatic {input.germline} {input.annotation} {input.combinedList} {input.hotspot} > {output.germline}
	else
		touch {input.germline}.dummy
		perl {input.convertor} germline {input.germline}.dummy {input.germline} {input.annotation} {input.combinedList} {input.hotspot} > {output.germline}.gl
		perl {input.convertor} somatic  {input.hotspot} {input.cgc} {input.combinedList}  {input.germline} {input.annotation} >{output.germline}.som
		perl {input.combine} {output.germline}.gl {output.germline}.som  >{output.germline}
		rm -rf {output.germline}.gl {output.germline}.som {input.germline}.dummy
	fi
	#######################
	"""
############
#       Actionable
############
rule Actionable_RNAseq:
	input:
		rnaseq    ="{subject}/{TIME,[0-9]+}/{subject}/db/{subject}.rnaseq",
		annotation="{subject}/annotation/{subject}.Annotations.coding.rare.txt",
		annotate  =NGS_PIPELINE + "/scripts/addAnnotations2vcf.pl",
		convertor =NGS_PIPELINE + "/scripts/" + config["Actionable_mutation"],
		cancerGeneCensus = config["annovar_data"]+"geneLists/CGCensus_Hereditary.txt",
		hotspot= config["annovar_data"]+"hg19_SomaticActionableSites.txt",
		combinedList  = config["annovar_data"]+"geneLists/combinedList_04292016",
		combine=NGS_PIPELINE + "/scripts/germlineOnly.pl",
		cgc           = config["annovar_data"]+"geneLists/CancerGeneCensus.v76.txt"
	output:
		rnaseq="{subject}/{TIME,[0-9]+}{ACT_DIR}{subject}.rnaseq.actionable.txt",
	params:
		rulename  = "ActionableMutations",
		batch    = config[config['host']]['job_default']
	shell: """
	#######################
	touch {input.rnaseq}.dummy
	perl {input.convertor} rnaseq {input.rnaseq}.dummy {input.rnaseq} {input.annotation} {input.combinedList} {input.hotspot} > {output.rnaseq}.gl
	perl {input.convertor} somatic  {input.hotspot} {input.cgc} {input.combinedList}  {input.rnaseq} {input.annotation} >{output.rnaseq}.som
	perl {input.combine} {output.rnaseq}.gl {output.rnaseq}.som  >{output.rnaseq}
	rm -rf {output.rnaseq}.gl {output.rnaseq}.som {input.rnaseq}.dummy
	#######################
	"""
############
#	**END**
############
