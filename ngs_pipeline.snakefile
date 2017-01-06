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
NOW=os.environ['TMP']
configfile: NGS_PIPELINE +"/config/config_annotation.json"
configfile: NGS_PIPELINE +"/config/config_common.json"
configfile: NGS_PIPELINE +"/config/config_cluster.json"
if HOST == 'biowulf.nih.gov':
	configfile: NGS_PIPELINE +"/config/config_common_biowulf.json"
elif HOST == 'login01':
	configfile: NGS_PIPELINE +"/config/config_common_tgen.json"

config['host'] = HOST
GROUP=config['group']
MAIL=config['mail']
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
sleep 20s
if [ {HOST} == 'biowulf.nih.gov' ]
	then
		MEM=`echo "${{SLURM_MEM_PER_NODE}} / 1024 "|bc`
		LOCAL="/lscratch/${{SLURM_JOBID}}/"
		THREADS=${{SLURM_CPUS_ON_NODE}}
elif [ {HOST} == 'login01' ]
	then
		module load slurm
		module load gcc/4.8.1
		MEM=`scontrol show job ${{SLURM_JOB_ID}} | grep "MinMemoryNode"| perl -n -e'/MinMemoryNode=(\d*)G/ && print $1'`
		mkdir -p /projects/scratch/ngs_pipeline_{NOW}_${{SLURM_JOB_ID}}/
		LOCAL="/projects/scratch/ngs_pipeline_{NOW}_${{SLURM_JOB_ID}}/"
		THREADS=`scontrol show job ${{SLURM_JOB_ID}} | grep  "MinCPUsNode" | perl -n -e'/MinCPUsNode=(\d*)/ && print $1'`
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
#make dictionary containing fastq file locations.
# Die if a library is ran twice or more.
FQ={}
for sample in config['library'].keys():
	for fq in config['library'][sample]:
		if len(config['library'][sample]) == 1:
			if os.path.isfile(DATA_DIR+fq+"/"+fq+"_R1.fastq.gz"):
				FQ[sample] =[DATA_DIR+fq+"/"+fq+"_R1.fastq.gz", DATA_DIR+fq+"/"+fq+"_R2.fastq.gz"]
			elif os.path.isfile(DATA_DIR+"Sample_"+fq+"/Sample_"+fq+"_R1.fastq.gz"):
				FQ[sample] =[DATA_DIR+"Sample_"+fq+"/Sample_"+fq+"_R1.fastq.gz", DATA_DIR+"Sample_"+fq+"/Sample_"+fq+"_R2.fastq.gz"]
			else:
				print("can not locate fastq file for sample", fq)
	  			exit()
		else:
			exit()
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
SUB_HOT = {}
SUB_IGV = {}
SUB_CON_QC = {}
SAMPLES =[]
somaticPairs = {}
somaticCopy = {}
SequenzaPairs ={}
pairedCapture = {}
# Inputs for the targets, where direct list can not be used.
for subject in config['subject'].keys():
	SUBS.append(subject)
	PATIENTS.append(subject)
	SUB_BAMS[subject]= ["{subject}/{TIME}/{sample}/{sample}.bwa.final.bam".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_COV[subject] = ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.coverage.txt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_HOT[subject] = ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.hotspot.depth".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_LOH[subject] = ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.loh".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_GT[subject]  = ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.gt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
	SUB_CON_QC[subject]  = ["{subject}/{TIME}/{sample}/qc/{sample}.consolidated_QC".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in config['subject'][subject]]
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
# Many of the targets.
ALL_FASTQC  = ["{subject}/{TIME}/{sample}/qc/fastqc/{sample}_R2_fastqc.html".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC      = ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.flagstat.txt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC     += ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.squeeze.done".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
#ALL_QC     += ["{subject}/{TIME}/{sample}/verifyBamID/{sample}.selfSM".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC     += ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.hotspot.depth".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
ALL_QC     += ["{subject}/{TIME}/{sample}/qc/{sample}.bwa.gt".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
#ALL_QC     += ["{subject}/{TIME}/{sample}/qc/BamQC/qualimapReport.html".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[s], sample=s) for s in SAMPLES]
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
ALL_QC     += expand("{subject}/{TIME}/annotation/AnnotationInput.coding.rare.txt", TIME=TIME, subject=PATIENTS)
ALL_QC     += expand("{subject}/{TIME}/annotation/{subject}.Annotations.coding.rare.txt", TIME=TIME, subject=PATIENTS)
ALL_QC     += expand("{subject}/{TIME}/qc/{subject}.config.txt", TIME=TIME, subject=PATIENTS)
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
			if config['sample_captures'][Tumor] not in config['Panel_List']:
				SequenzaPairs[Tumor] = ["{subject}/{TIME}/{sample}/{sample}.mpileup.gz".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[Normal], sample=Normal), "{subject}/{TIME}/{sample}/{sample}.mpileup.gz".format(TIME=TIME, subject=SAMPLE_TO_SUBJECT[Tumor], sample=Tumor) ]
###########################################################################
# This is to make list of DB file list. (germline, variants, somatic, rnaseq)
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
###########################################################################
# This is to find out if we need to make variants db file or germline file.
ACT_TYPE =[]
DECIDE_GL={}
for subject in config['subject'].keys():
	normal = None
	tumor  = None
	pair   = None
	for sample in config['subject'][subject]:
		if config['sample_type'][sample] == 'Tumor':
			tumor = 'yes'
			if sample in config['sample_references'].keys():
				pair  = 'yes'
		elif config['sample_type'][sample] == 'Normal':
			normal = 'yes'
	if pair =='yes':
		DECIDE_GL[subject] = 'gl_only'
	elif pair == None:
		if tumor == None and normal =='yes':
			DECIDE_GL[subject] = 'gl_only'
		else:
			ACT_TYPE +=[subject]
###########################################################################
# To make union somatic file for every library. 
#	This includes logic to get files based on sample_ref and sample_rnaseq
UNION_SOM_MUT={}
UNION_SOM_MUT_LIST =[]
for sample in config['sample_references'].keys():
	subject=SAMPLE_TO_SUBJECT[sample]
	local =[(subject+"/"+TIME+"/"+sample+"/calls/"+sample+".strelka.snvs.annotated.txt"),
		(subject+"/"+TIME+"/"+sample+"/calls/"+sample+".strelka.indels.annotated.txt"),
		(subject+"/"+TIME+"/"+sample+"/calls/"+sample+".MuTect.annotated.txt")]
	if sample in config['sample_RNASeq'].keys():
		local = [w.replace('annotated','annotated.expressed') for w in local]
	UNION_SOM_MUT[sample] = local
	UNION_SOM_MUT_LIST +=[subject+"/"+TIME+ACT_DIR+sample+".unionSomaticVars.txt"]
	if config['sample_captures'][sample] not in config['Panel_List']:
		UNION_SOM_MUT_LIST +=[subject+"/"+TIME+ACT_DIR+sample+".mutationalSignature.pdf"]

##########################################################################
# To create lists to be filled in SUBJECT_ANNO
for subject in config['subject']:
	local  = []
	for sample in config['subject'][subject]:
		local.extend([(subject+"/"+TIME+"/"+sample+"/calls/"+sample+".HC_DNASeq.snpEff.txt"),
			      (subject+"/"+TIME+"/"+sample+"/calls/"+sample+".Platypus.snpEff.txt"),
			      (subject+"/"+TIME+"/"+sample+"/calls/"+sample+".bam2mpg.snpEff.txt")])
	if subject not in SUBJECT_VCFS:
		SUBJECT_VCFS[subject] = local
	if subject in ACT_TYPE:
		germline = [w.replace('snpEff','annotated') for w in local]
		add_to_SUBJECT_ANNO(subject,"variants",germline)
	else:
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
	if config['sample_captures'][sample] not in config['Panel_List']:
		COPY_NUMBER +=[subject+"/"+TIME+"/"+sample+"/sequenza/"+sample+"/"+sample+"_alternative_fit.pdf"]
		COPY_NUMBER +=[subject+"/"+TIME+"/"+sample+"/sequenza/"+sample+".txt"]
		COPY_NUMBER +=[subject+"/"+TIME+"/"+sample+"/NeoAntigen/"+sample+".somatic.vep.vcf"]
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
# Expressed Mutations
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
localrules: Khanlab_Pipeline, RNASeq
#IGV_Session, DBinput, AttachAnnotation, Expressed, vcf2txt, symlink_tophatBam, copyNovoBam, Actionable_Germline, Actionable_RNAseq, Actionable_Somatic, Actionable_Variants, Actionable_fusion, Sub_Fusion, makeConfig, TargetInterval, QC_Summary_Patient,QC_Summary,UnionSomaticCalls,TOPHAT_LINK, SampleGT,QC_Sum, FormatInput, RNASeqQC_1,RNASeqQC1 RNASeqQC_2,RNASeqQC_3, Cuff_Mat
#Circos, CoveragePlot, BoxPlot_Hotspot, makeConfig,Ideogram
###########################################################################
#                               Rule Book				  #
###########################################################################
include: NGS_PIPELINE +"/ruleBook/bamUtil.snakefile"
include: NGS_PIPELINE +"/ruleBook/verifyBamID.snakefile"
include: NGS_PIPELINE +"/ruleBook/rnaseq_pipeline.snakefile"
include: NGS_PIPELINE +"/ruleBook/rnaseqQC.snakefile"
include: NGS_PIPELINE +"/ruleBook/readDepth.snakefile"
include: NGS_PIPELINE +"/ruleBook/failedExon.snakefile"
include: NGS_PIPELINE +"/ruleBook/hsMetrix.snakefile"
include: NGS_PIPELINE +"/ruleBook/Consolidate.snakefile"
include: NGS_PIPELINE +"/ruleBook/universal.snakefile"
include: NGS_PIPELINE +"/ruleBook/mutationalSignature.snakefile"
include: NGS_PIPELINE +"/ruleBook/NeoAntigen.snakefile"

include: NGS_PIPELINE +"/ruleBook/haplotypeCaller.snakefile"
include: NGS_PIPELINE +"/ruleBook/platypus.snakefile"
include: NGS_PIPELINE +"/ruleBook/bam2mpg.snakefile"

include: NGS_PIPELINE +"/ruleBook/gatk_RNASeq.snakefile"
include: NGS_PIPELINE +"/ruleBook/ideogram.snakefile"
include: NGS_PIPELINE +"/ruleBook/Actionable.snakefile"
include: NGS_PIPELINE +"/ruleBook/UnionSomaticMutations.snakefile"
include: NGS_PIPELINE +"/ruleBook/plots.snakefile"
include: NGS_PIPELINE +"/ruleBook/annot.snakefile"
include: NGS_PIPELINE +"/ruleBook/STAR.snakefile"

include: NGS_PIPELINE +"/ruleBook/Sequenza.snakefile"


ALL_VCFs =[]
for subject in SUBJECT_VCFS.keys():
	for vcf in SUBJECT_VCFS[subject]:
		vcf = vcf.replace('snpEff.txt', 'raw.vcf')
		ALL_VCFs +=[vcf]
		vcf = vcf.replace('raw.vcf', 'raw.snpEff.vcf')
		ALL_VCFs +=[vcf]
###########################################################################
onerror:
	shell("find .snakemake/ ! -readable -prune \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
        shell("find .snakemake/ ! -readable -prune -group $USER -exec chgrp -f {GROUP} {{}} \;")
	shell("find {PATIENTS} -group $USER -exec chgrp -f {GROUP} {{}} \;")
	shell("find {PATIENTS} \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
	shell("ssh {HOST} \"echo 'Pipeline failed on {PATIENTS}. Error occured on {HOST}. Working Dir:  {WORK_DIR}' |mutt -s 'Khanlab ngs-pipeline Status' `whoami`@mail.nih.gov  {MAIL} \"")
	shell("find .snakemake/ ! -readable -prune -group $USER -exec chgrp -f {GROUP} {{}} \;")
onstart:
	f = open('ngs_pipeline_%s.csv' % NOW , 'w')
	print ('#Patient','Diagnosis','CaseID',sep='\t', end='\n',file=f)
	for subject in sorted(PATIENTS):
		diagnosis =config['Diagnosis'][SUBJECT_TO_SAMPLE[subject][0]]
		print (subject,diagnosis,TIME,sep='\t', end='\n',file=f)
	
	shell("for sub in {PATIENTS}; do rm -rf {WORK_DIR}/${{sub}}/{TIME}/successful.txt ; done")
	shell("ssh {HOST} \"echo 'ngs-pipeline started on {PATIENTS} on {HOST}. Working Dir:  {WORK_DIR}' |mutt -s 'Khanlab ngs-pipeline Status' `whoami`@mail.nih.gov {MAIL} \"")
onsuccess:
	shell("find .snakemake/ ! -readable -prune \( -type f -user $USER -exec chmod g+r {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)")
	shell("find .snakemake/ ! -readable -prune -group $USER -exec chgrp -f {GROUP} {{}} \;")
	print("Workflow finished, no error")
###########################################################################
rule Khanlab_Pipeline:
	input:
		SUB_IGV.values(),
		COPY_NUMBER,
		ALL_VCFs,
		CON_QC,
		ALL_QC,
		ALL_FASTQC,
		varFiles,
		DBFiles,
		ActionableFiles,
		UNION_SOM_MUT_LIST
	version: "1.0"
	wildcard_constraints:
		NOW="\w+"
	params:
		rulename = "Final",
		batch    = config[config['host']]["job_default"],
		group    = config["group"],
		wait4job = NGS_PIPELINE + "/scripts/block_for_jobid.pl",
		sort 	 = NGS_PIPELINE + "/scripts/awk_sort_withHeader.awk",
		mail 	 = NGS_PIPELINE + "/scripts/tsv2html.final.sh",
		email    = config["mail"],
		host     = config["host"],
		subs     = PATIENTS
	shell: """
	#######################
	find {PATIENTS} log -group $USER -exec chgrp -f {params.group} {{}} \;
	find {PATIENTS} log \( -type f -user $USER -exec chmod g+rw {{}} \; \) , \( -type d -user $USER -exec chmod g+rwx {{}} \; \)
	export LC_ALL=C
	
	for sub in {params.subs}
        do
                touch {WORK_DIR}/${{sub}}/{TIME}/successful.txt
		chmod g+rw {WORK_DIR}/${{sub}}/{TIME}/successful.txt 	
		chgrp {params.group} {WORK_DIR}/${{sub}}/{TIME}/successful.txt
        done
	ssh {params.host} "{params.mail} --location {WORK_DIR} --host {params.host} --head {WORK_DIR}/ngs_pipeline_{NOW}.csv |mutt -e \\\"my_hdr Content-Type: text/html\\\" -s 'Khanlab ngs-pipeline Status' `whoami`@mail.nih.gov {params.email}"
	rm -rf {WORK_DIR}/ngs_pipeline_{NOW}.csv
	#######################
	"""
############
# Print Config to a file
############
rule makeConfig:
	output:	"{subject}/{TIME}/qc/{subject}.config.txt" 
	params:
		rulename = "configPrint",
		batch    = config[config['host']]["job_default"],
		hash     = json.dumps(config, sort_keys=True)
	shell: """
	#######################
	echo '{params.hash}'  >{output}
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
	ssh {params.host} "sh {input.mail} --name {wildcards.subject} --diagnosis '{params.diagnosis}' --head {WORK_DIR}/{wildcards.subject}/{TIME}/qc/{wildcards.subject}.genotyping.txt | mutt -e \\\"my_hdr Content-Type: text/html\\\" -s 'Genotyping Result on {wildcards.subject}' `whoami`@mail.nih.gov {params.mail} "
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
#       Somatic Copy Number LRR (Median Corrected)
############
rule CN_LRR:
	input:
		files=lambda wildcards: somaticCopy[wildcards.Tumor],
		ref=config["gene_coord"],
		index=config["reference"].replace('.fasta', '.index.txt'),
		tool=NGS_PIPELINE+ "/scripts/AddGene.pl",
		cgc    = config["annovar_data"]+config["geneList"],
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
		interval=lambda wildcards: config['hotspot_bed'][config['sample_captures'][wildcards.sample]],
		genome=config["reference"].replace(".fasta",".genome")
	output: "{base}/{TIME}/{sample}/qc/{sample}.{aligner}.hotspot.depth"
	version: config["bedtools"]
	params:
		rulename  = "HotSpotCov",
		batch     = config[config['host']]["job_hotspot"],
	shell: """
	#######################
	module load samtools
	module load bedtools/{version}
	slopBed -i {input.interval} -g {input.genome} -b 50 >${{LOCAL}}/Region.bed
	samtools view -hF 0x400 -q 30 -L ${{LOCAL}}/Region.bed {input.bam} | samtools view -ShF 0x4 - | samtools view -SuF 0x200 - | bedtools coverage -abam - -b {input.interval} >{output}
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
	dir=`echo {params.work_dir} | sed -e 's/\/data\/khanlab/K:/g' | sed -e 's/\/projects\/Clinomics/Y:/g' |sed -e 's/\/data\/Clinomics/V:/g'`
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
	#######################
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
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T IndelRealigner -R {input.ref} -known {input.phase1} -known {input.mills} -I {input.bam} --targetIntervals ${{LOCAL}}/{wildcards.sample}.realignment.intervals -o ${{LOCAL}}/{wildcards.sample}.lr.bam --maxReadsInMemory 1500000
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T BaseRecalibrator -R {input.ref} -knownSites {input.phase1} -knownSites {input.mills} -I ${{LOCAL}}/{wildcards.sample}.lr.bam -o ${{LOCAL}}/{wildcards.sample}.recalibration.matrix.txt
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T PrintReads -R {input.ref} -I ${{LOCAL}}/{wildcards.sample}.lr.bam -o {output.bam} -BQSR ${{LOCAL}}/{wildcards.sample}.recalibration.matrix.txt
	mv -f {wildcards.base}/{TIME}/{wildcards.sample}/{wildcards.sample}.bwa.final.bai {output.index}
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
		interval=lambda wildcards: config['target_intervals'][pairedCapture[wildcards.Tumor]].replace("target","targetbp")
	output:
		snps="{subject}/{TIME}/{Tumor}/calls/{Tumor}.strelka.snvs.raw.vcf",
		indels="{subject}/{TIME}/{Tumor}/calls/{Tumor}.strelka.indels.raw.vcf"
	version: config["strelka"]
	params:
		rulename = "Strelka",
		batch    = config[config['host']]["job_strelka"],
		config=NGS_PIPELINE + "/Tools_config/"+config["strelka_config"],
		vcftools = config["vcftools"]
	shell: """
	#######################
	module load strelka/{version}
	configureStrelkaWorkflow.pl --normal={input[2]} --tumor={input[0]}\
	--ref={input.ref} --config={params.config} --output-dir=${{LOCAL}}/strelka
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
#	snpEff
############
rule SNPEff:
	input:
		vcf="{subject}/{TIME}/{sample}/calls/{base}.raw.vcf",
		ref=config["reference"],
	output:
		eff="{subject}/{TIME}/{sample}/calls/{base}.raw.snpEff.vcf"
	version: config["snpEff"]
	params:
		rulename      ="snpEff",
		batch	      =config[config['host']]["job_snpeff"],
		snpEff_genome =config["snpEff_genome"],
		snpEff_config=NGS_PIPELINE + "/Tools_config/"+config["snpEff_config"],
		annovar       =config["annovar"]
	shell: """
	#######################
	module load snpEff/{version}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $SNPEFF_JARPATH/SnpSift.jar dbnsfp -c {params.snpEff_config} -a {input.vcf} | java -Xmx${{MEM}}g -jar $SNPEFF_JARPATH/snpEff.jar -t -canon {params.snpEff_genome} > {output.eff}
	#######################
	"""
############
#       vcf2txt
############
rule VCF2TXT:
	input:
		eff="{subject}/{TIME}/{sample}/calls/{base}.raw.snpEff.vcf",
		vcf2txt=NGS_PIPELINE + "/scripts/vcf2txt.pl"
	output:
		txt="{subject}/{TIME}/{sample}/calls/{base}.snpEff.txt"
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
		temp("{subject}/{TIME}/annotation/AnnotationInput.anno"),
		temp("{subject}/{TIME}/annotation/AnnotationInput.sift")
	version: config["annovar"]
	params:
		rulename   = "FormatInput",
		batch      = config[config['host']]["job_default"],
		fAEV       = NGS_PIPELINE + "/scripts/findAlreadyExistingVariants.pl"
	shell: """
	#######################
	module load annovar/{version}
	cut -f 1-5 {input.txtFiles} |sort |uniq > {wildcards.subject}/{TIME}/annotation/AnnotationInput
	perl {input.convertor} {wildcards.subject}/{TIME}/annotation/AnnotationInput
	rm -rf "{wildcards.subject}/{TIME}/annotation/AnnotationInput.pph"
	#######################
	"""
############
#       CopyAnnotationFile
############
rule CopyAnnotationFile:
	input:
		"{subject}/{TIME}/annotation/AnnotationInput.coding.rare.txt"
	output:
		"{subject}/{TIME}/annotation/{subject}.Annotations.coding.rare.txt"
	params:
		rulename   = "caf",
		batch      = config[config['host']]["job_default"],
	shell: """
	#######################
	cp {input} {output}
	#######################
	"""
############
#	Add Annotation back to sample level file
############
rule AttachAnnotation:
	input:
		txt="{subject}/{TIME}/{base1}/calls/{base}.snpEff.txt",
		ref="{subject}/{TIME}/annotation/{subject}.Annotations.coding.rare.txt",
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
		Mutation="{subject}/{TIME}/{sample}/calls/{base}.annotated.txt",
		convertor = NGS_PIPELINE + "/scripts/mpileup.pl"
	output: "{subject}/{TIME}/{sample}/calls/{base}.annotated.expressed.txt"
	version: config["samtools"]
	params:
		rulename  = "Expressed",
		batch     = config[config['host']]["job_expressed"],
		name      = lambda wildcards: config["sample_RNASeq"][wildcards.sample]
	shell: """
	#######################
	module load samtools/{version}
	perl {input.convertor} {input.Mutation} {wildcards.subject}/{TIME}/{params.name}/{params.name}.star.final.bam {input.RNASeq} >{output}	
	#######################
	"""
############
#       Database Input
############
rule DBinput:
	input:
		txtFiles=lambda wildcards: SUBJECT_ANNO[wildcards.subject][wildcards.group],
		convertor=NGS_PIPELINE + "/scripts/makeDBVariantFile.pl",
		tool=NGS_PIPELINE + "/scripts/AddSampleType.pl",
		tool1=NGS_PIPELINE + "/scripts/addFS.pl",
		txtFiles1=lambda wildcards: SUBJECT_VCFS[wildcards.subject]
	output: "{subject}/{TIME}/{subject}/db/{subject}.{group}"
	params:
		rulename = "makeDBinput",
		batch    = config[config['host']]['job_default'],
		hash 	 = lambda wc: " ".join("{} {}".format(a, b) for a, b in config["sample_type"].items()),
		hash1	 = lambda wc: " ".join("{} {}".format(a, b) for a, b in config["sample_captures"].items()),
	shell: """
	#######################
	perl {input.convertor} {input.txtFiles} |perl {input.tool} - "{params.hash}" "{params.hash1}" >{output}.tmp
	if [ {wildcards.group}  == 'germline' ]; then	
		perl {input.tool1} {output}.tmp {input.txtFiles1} >{output}
		rm -rf {output}.tmp
	else
		mv {output}.tmp {output}
	fi
	#######################
	"""
