TARGET     += expand("{subject}/{TIME}/qc/{subject}.coveragePlot.png",TIME=TIME, subject=PATIENTS)
TARGET     += expand("{subject}/{TIME}/qc/{subject}.hotspot_coverage.png", TIME=TIME, subject=PATIENTS)
if 'subject' in config:
	for subject in config['subject']:
		for library in config['subject'][subject]:
			if config['sample_captures'][library] not in config['Panel_List']:
				# any output which is desired on all libraries but Panel goes here, the list of panel captures should be maintained in the Panel_List in config file
				TARGET    += [subject+"/"+TIME+"/qc/"+subject+".circos.png"]
############
# CoveragePlot
############
rule CoveragePlot:
	input:
		covFiles=lambda wildcards: SUB_COV[wildcards.subject],
		coverage =NGS_PIPELINE + "/scripts/coverage.R"
	output: "{subject}/{TIME}/qc/{subject}.coveragePlot.png",
	version: config["version_R"]
	params:
		rulename = "covplot",
		batch    = config[config['host']]["job_covplot"]
	shell: """
	#######################

	cp -f {input.covFiles} ${{LOCAL}}

	module load R/{version}
	R --vanilla --slave --silent --args ${{LOCAL}} {output} {wildcards.subject} <{input.coverage}
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
	version: config["version_R"]
	params:
		rulename = "Circos",
		batch    = config[config['host']]["job_covplot"]
	shell: """
	#######################
	cp -f {input.lohFiles} ${{LOCAL}}
	module load R/{version}
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
	version: config["version_R"]
	params:
		rulename = "Boxplot",
		batch    = config[config['host']]["job_covplot"]
	shell: """
	#######################
	cp -f {input.covFiles} ${{LOCAL}}
	module load R/{version}
	R --vanilla --slave --silent --args ${{LOCAL}} {output} {wildcards.subject} <{input.boxplot}
	#######################
	"""
