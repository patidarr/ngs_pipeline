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
	module load vcftools/{params.vcftools} {params.samtools}
	for CHR in `seq 1 22` X Y;
	do
	echo "bam2mpg --shorten --vcf_spec --qual_filter 20 -bam_filter '-q31' --region chr${{CHR}} --only_nonref --snv_vcf ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.snps.vcf --div_vcf ${{LOCAL}}/chr${{CHR}}{wildcards.sample}.indel.vcf {input.ref} {input.bam}"
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
