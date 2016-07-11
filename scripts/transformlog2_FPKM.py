#! /usr/bin/python

import sys 
import math
from math import sqrt,log
import re

header=''

def make_dictionary(file_array):

	cufflinks_dict={}

        for i in range(0,len(file_array)):
		if file_array[i] != '':
			prelim_info_list=[]
                	each_gene_list=file_array[i].split("\t")
						
			try:	
				##Prelimnary Info
				entry_name=each_gene_list[0]
				Gene_ID=each_gene_list[3]
				Gene_Name=each_gene_list[4]
				tss_id=each_gene_list[5]
				locus=each_gene_list[6]
				length=each_gene_list[7]
				coverage=each_gene_list[8]
				FPKM=each_gene_list[9]
				
				log2_FPKM=math.log(float(FPKM)+1)/math.log(2)
				
				if entry_name not in cufflinks_dict:
					cufflinks_dict[entry_name]=entry_name+"\t"+Gene_ID+"\t"+Gene_Name+"\t"+tss_id+"\t"+locus+"\t"+str("%.4f" % log2_FPKM)
				else:
					pass
			except:
				pass
	
	return cufflinks_dict
def print_dict(exp_type,exon_dict_2):
	
	global header
	main_string=''
	key_set=exon_dict_2.keys()
	
	if exp_type=="genes":
		header ="Tracking_ID"+"\t"+"Gene_id"+"\t"+"Gene_Symbol"+"\t"+"tss_id"+"\t"+"locus"+"\t"+"log2(FPKM)"+"\n"
	elif exp_type=="isoforms":
		header="Gene_id"+"\t"+"Gene_Symbol"+"\t"+"Transcript_id"+"\t"+"tss_id"+"\t"+"locus"+"\t"+"log2(FPKM)"+"\n"
	
	for key in key_set:
		exon_details=exon_dict_2[key]
		main_string += exon_details+"\n"
	
	return main_string
	
if __name__ == '__main__' :

	output_dict={}
	
	expression_type=sys.argv[1]
	input_file=sys.argv[2]
	output_file=sys.argv[3]
        
	ref_file=open(input_file,"r")
	ref_file_data=ref_file.read()
        ref_data=ref_file_data.split("\n")
        output_dict=make_dictionary(ref_data)

	to_print=print_dict(expression_type,output_dict)
		
	save_results=open(output_file,'wb')
        save_results.write(header) 	
	save_results.write(to_print)
