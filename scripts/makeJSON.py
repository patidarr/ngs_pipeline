#!/usr/bin/python
import os
import json
from sys import argv

f = open(argv[1], 'r')
samples=argv[2]
output ={}
output["sample_captures"] ={}
output["Diagnosis"] = {}
output["sample_type"] = {}
output["subject"] = {}
output["library"] = {}
output["sample_references"] = {}
output["sample_RNASeq"] = {}
output["RNASeq"] = {}

patientIndex=0
TypeIndex=1
DiagnosisIndex=2
captureIndex=4
FCIDIndex=5
libraryIndex=6
normRefIndex=7
rnaRefIndex=8
#ID
#custom ID
#Type
#Diagnosis
#Type of sequencing
#Enrichment step
#FCID
#Library ID
#Matched normal
#Matched RNA-seq lib
#Case Name

for line in f:
	line = line.rstrip()	
	column = line.split("\t")
	column[captureIndex] =column[captureIndex].lower()
	for sample in samples.split(','):
		if column[patientIndex] == sample:
			output["sample_captures"][column[libraryIndex]]=column[captureIndex]
			output["Diagnosis"][column[libraryIndex]]=column[DiagnosisIndex]
			output["sample_type"][column[libraryIndex]]=column[TypeIndex]
			if column[normRefIndex]:
				output["sample_references"][column[libraryIndex]]=[column[normRefIndex]]
			if column[rnaRefIndex]:
				output["sample_RNASeq"][column[libraryIndex]] =[column[rnaRefIndex]]
			if not column[FCIDIndex]:
				output["library"][column[libraryIndex]] = [column[libraryIndex]]
			else:
				output["library"][column[libraryIndex]] = [column[libraryIndex]+"_"+column[FCIDIndex]]
			if 'DNA' in column[TypeIndex]:
				if column[patientIndex] not in output["subject"].keys():	
					output["subject"][column[patientIndex]]=[column[libraryIndex]]
				else:
					output["subject"][column[patientIndex]].append(column[libraryIndex])
			else:
				if column[patientIndex] not in output["RNASeq"].keys():
					output["RNASeq"][column[patientIndex]]=[column[libraryIndex]]
				else:
					output["RNASeq"][column[patientIndex]].append(column[libraryIndex])
		

print(json.dumps(output, sort_keys=True, indent=4))
