############
#       Union Somatic Calls
############
rule UnionSomaticCalls:
        input:
		txtFiles=lambda wildcards: SUBJECT_ANNO[wildcards.subject][wildcards.group],
		tool=NGS_PIPELINE + "/scripts/UnionSomaticCalls.pl"
        output: "{subject}/{TIME,[0-9]+}/{ACT_DIR}{sample}.uniqueSomaticCalls.txt"
        params:
                rulename = "UnionSomaticCalls",
                batch    = config[config['host']]['job_default']
        shell: """
        #######################
        perl {input.tool} {input.strelkaSNP} {input.strelkaIndel} {input.Mutect} >{output}
        #######################
        """

