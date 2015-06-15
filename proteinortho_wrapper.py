#!/usr/bin/python

#########################################################################################################################
## This is a python wrapper script for the tool Proteinortho (Lechner et al., BMC Bioinformatics, 2011). The script here 
## is used with specific parameters "-selfblast -identity=80 -cpus=8 -cov=80". 
########################################################################################################################

import argparse
import os
import sys
import subprocess
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import random

# use argparse to run through the command line options given
args_parser = argparse.ArgumentParser(description="Program for creating phylogenes the concatenated amino acid sequences of marker genes", epilog="Center for Microbial Oceanography Research and Education")
args_parser.add_argument('-i', '--input_folder', required=True, help='Folder of protein files with .faa at the end.')
#args_parser.add_argument('-o', '--output_folder', required=True, help='output folder.')
#args_parser.add_argument('-d', '--db', required=True, help='Database for matching.')
args_parser = args_parser.parse_args()

# set up object names for input/output/database folders
input_dir = args_parser.input_folder

#####################################
# Make directories to deposit results
#####################################
dirs = os.listdir(input_dir)
if "ortholog_clusters" in dirs:
	pass
else:
	cmd = "mkdir " + input_dir + "/ortholog_clusters"
	subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
	cmd = "mkdir " + input_dir + "/nr_proteins"
	subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
	cmd = "mkdir " + input_dir + "/nr_proteins/results"
	subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))

##################################################################################################
# Loop through protein files and make sure they are non-redundant
##################################################################################################
protein_dict = {}
protein_tally = {}
record_dict = {}
dirs = os.listdir(input_dir)
for filenames in dirs:
	if filenames.endswith(".faa"):
		handle = open(input_dir+'/'+filenames, "rU")
		o = open(input_dir+'/nr_proteins/'+filenames+".nr", 'w')

		for proteins in SeqIO.parse(handle, "fasta"):
		
		#record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
		#handle.close()

		#o = open(input_dir+filenames+".nr", 'w')
		#for proteins in record_dict:
			prot_id = proteins.id
			prot_seq = proteins.seq
			prot_desc = proteins.description
			list2 = prot_desc.split(' ')
			prot_desc = " ".join(list2[1:len(list2)])
			prot_length = len(prot_seq)

			if prot_id in protein_tally:
				protein_tally[prot_id] = protein_tally[prot_id] + 1
				new_id = prot_id+'_'+str(protein_tally[prot_id])
				record = SeqRecord(seq=prot_seq, id=new_id, description=prot_desc)
				protein_dict[new_id] = record
			else:
				protein_tally[prot_id] = 1
				new_id = prot_id+'_'+str(protein_tally[prot_id])
				record = SeqRecord(seq=prot_seq, id=new_id, description=prot_desc)
				protein_dict[new_id] = record
			SeqIO.write(protein_dict[new_id], o, "fasta")
			#o.write(str(prot_id)+' '+str(prot_desc)+' '+str(prot_length)+'\n')
			#o.write(">"+new_id+' '+prot_desc+'\n'+str(prot_seq)+'\n')
		o.close

##################################################################################################
# Run proteinortho on non-redunant protein files
##################################################################################################

#file_list = []
#directory = '/slipstream/home/faylward/python/proteinortho_wrapper/'+input_dir+'/nr_proteins'
#directory = '/slipstream/home/faylward/python/proteinortho_wrapper/'+input_dir+'/nr_proteins'
file_list = str()
dirs = os.listdir(input_dir+'/nr_proteins')
for files in dirs:
	if files.endswith(".nr"):
		#file_list = str(file_list)+' '+input_dir+'/nr_proteins/'+str(files)
		file_list = str(file_list)+' '+str(files)
		cmd = "formatdb -p T -i " + files
		print cmd
		#subprocess.call(cmd, shell=True, cwd=directory)
		subprocess.call(cmd, shell=True, cwd=input_dir+"/nr_proteins")

		#file_list.append(input_dir+'/nr_proteins/'+files)
#print file_list
cmd = "proteinortho -p=blastp -project="+input_dir+" -selfblast -identity=80 -cpus=8 -cov=80 "+file_list
print cmd
#directory = $HOME/python/proteinortho_wrapper/
#directory = '/slipstream/home/faylward/python/proteinortho_wrapper/test_proteins/nr_proteins'
subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'wb'), cwd=input_dir+"/nr_proteins")

cmd = "proteinortho5_singletons " + file_list + " < " + input_dir+'.proteinortho >'+'results/'+input_dir+'.singletons'
print cmd
subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'wb'), cwd=input_dir+"/nr_proteins")

cmd = "cat "+input_dir+'.proteinortho '+'results/'+input_dir+'.singletons >'+'results/'+input_dir+'.proteinortho.final'
print cmd
subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'wb'), cwd=input_dir+"/nr_proteins")


#############################################################################################################################################################
# Sort through output and add cluster names and annotations, provide annotation for the longest protein in a given cluster (randomly chosen if more than one)
#############################################################################################################################################################
f = open(input_dir+'/nr_proteins/results/'+input_dir+".proteinortho.final", 'r')
o = open(input_dir+'/nr_proteins/results/'+input_dir+".proteinortho.final.annote", 'w')
proteins = []
i = 0
for line in f.readlines():
	if line.startswith('#'):
		o.write("Cluster\tAnnotation\t"+line)
	else:
		proteins = []
		i +=1
		list1 = line.split('\t')
		list2 = list1[3:len(list1)]
		for items in list2:
			if '*' in items:
				pass
			elif ',' in items:
				list3 = items.split(',')
				for new_items in list3:
					proteins.append(new_items.rstrip())
			else:
				proteins.append(items.rstrip())
	prot_length = {}
	for prot in proteins:
		length = int(len(protein_dict[prot].seq))
		prot_length[prot] = length
		rand = random.choice([key for key,val in prot_length.iteritems() if val == max(prot_length.values())])
		description = protein_dict[rand].description
		cluster_id = input_dir+"_Cluster_"+str(i)
		o.write(cluster_id+'\t'+description+'\t'+line)
o.close()
		
		
	





