#!/usr/bin/python
import os
import json
from sys import argv

print("RuleName\tMinTime\tMaxTime\tMeanTime")
with open(argv[1]) as data_file:
	data =json.load(data_file)

for RULE in sorted(data['rules'].keys()):
	a = None
	b = None
	c = None
	for TIME in data['rules'][RULE]:
		if TIME == 'min-runtime':
			a=round((data['rules'][RULE][TIME]/60))
		elif TIME == 'max-runtime':
			b=round((data['rules'][RULE][TIME]/60))
		elif TIME == 'mean-runtime':
			c=round((data['rules'][RULE][TIME]/60))
	if c >0:
		print RULE,"\t",a,"\t",b,"\t",c
