############
##	Tiled data file(.tdf)
############
rule BAM2TDF:
	input:
		bam="{base}/{TIME}/{sample}/{sample}.{aligner}.final.bam",
		bam_bai="{base}/{TIME}/{sample}/{sample}.{aligner}.final.bam.bai"
	output:
		"{base}/{TIME}/{sample}/{sample}.{aligner}.final.bam.tdf"
	version: config['igvtools']
	params:
		rulename        = "bamtdf",
		ref             = config["reference"],
		batch           = config[config['host']]['job_igvtools']
	shell:  """
	#######################
	module load igvtools/{version}

	cd ${{LOCAL}}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} $IGVTOOLSHOME/igvtools.jar count {WORK_DIR}/{input.bam} {WORK_DIR}/{output} {params.ref}
	#######################
	"""
