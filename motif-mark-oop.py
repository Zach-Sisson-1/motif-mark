#!/usr/bin/env python

import argparse
from collections import defaultdict
from distutils.command.install_lib import PYTHON_SOURCE_EXTENSION
from operator import index
import re
import itertools
from urllib.request import parse_http_list

def get_args():
	parser = argparse.ArgumentParser(description="XXXX")
	parser.add_argument("-f", "--fasta_file", help = "Path to the input Fasta File", type=str, required=True)
	parser.add_argument("-m", "--motifs_file", help = "Path to the input Motifs file", type=str, required=True)
	return parser.parse_args()

args = get_args()


#Initalize class for Motifs:

class Uniq_Motif:
	"""
	Motif class that takes in a motif and can search for that motifs in a supplied gene, outputing coordinates
	"""
	def __init__(self, motif):
		self.motif = motif  #store the unaltered motifs
		self.motif_lower = motif.lower()   #generates list of lowercase motifs
		#Read each nucleotide in the motif and append the appropriate regex search rule to create unique regex
		regex = '(?=('
		for nucleotide in self.motif_lower:
			if nucleotide == "r":
				regex += "(a|g)"   						#considers purine a,,g
			elif nucleotide =="y":
				regex += "(c|t|u)"			 				#considers pyrimadine c,t,u
			elif nucleotide == "u":
				regex += "(t|u)"
			elif nucleotide == "t":
				regex += "(t|u)"
			else:										#else add whatever the nucelotide is to the running regex
				regex += nucleotide
		regex += '))'					
		self.regex = regex								#Assign unique regex expression
		print(regex)
	
	def gene_search(self,Gene):
		"""
		Takes in a gene object and searches for the motif in the gene, returning a dictionary  
		 containing the locations of the matching locations as values.
		""" 
		#initialize dict which will house motif as keys and location hits as values 
		match_dict = {}

		#Create iterobjects for the motif and locate all the motif hits within the given sequence 				
		match_dict[self.motif] = re.finditer(self.regex,Gene.gene_seq_lower)

		#Use only the the index locations of the matches
		index_dict = {}	
		for motif,iterobj in match_dict.items():
			lst = [obj.span() for obj in list(iterobj)]
			index_dict[motif] = lst
		
		return index_dict
	

class Gene:
	"""
	Gene class that takes in a gene pattern and stores all exons, pre-exons and post-exons.
	"""
	def __init__(self,gene_seq, header):
		self.gene_seq = gene_seq
		self.gene_seq_lower = gene_seq.lower()   #hold lowercase string for motif searching
		self.header = header						#header for FASTA seq 
		self.pre_exon = ""  
		self.exon = ""
		self.post_exon = ""	

	#Parse gene_seq into pre-exon, exon, and postexon
		self.pre_exon = re.findall('[a-z]+',self.gene_seq)[0]
		self.post_exon = re.findall('[a-z]+',self.gene_seq)[1]
		self.exon = re.findall('[A-Z]+',self.gene_seq)



#read in motif file, creating a list to hold motifs
motif_list = [] 
with open(args.motifs_file, "r") as motif_file:
	for line in motif_file:
		Motif = line.strip('\n')
		motif_list.append(str(Motif))




#Read in fasta file, creating a dict to hold gene objects
with open(args.fasta_file, "r") as fasta_file:
	gene_dict = {} #empty dict that will contain fasta headers and corresponding genes
	for line in fasta_file:
		if line[0] == '>':
			gene_name = re.search('[A-Z].+',line)[0]    #grabs the gene name from Fasta header"
			gene_dict[gene_name] = ""			
		else:
			gene_dict[gene_name] += line.strip('\n')			#builds string of entire gene sequence


# #######Test for creating motif list from file WORKS
# print("making sure motif file is being read correctly")
# print('\n')
# print("these are the detected motifs")
# print(motif_list)
# print('\n')



# ######Test for gene_dict         WORKS
# print("Testing Gene+seq Dictionary from Fasta File")
# for gene,sequence in gene_dict.items():
# 	print(gene)
# 	print('\n')
# 	print(sequence)
# 	print('\n')


#Create dict of Gene object for each gene in dictionary.   Keys = header name and values = gene objects 
Genes = {name:Gene(gene_seq=(gene_dict[name]), header=name) for name in gene_dict.keys()}

	
#test  for Gene objects   WORKS
#print(Genes)

###########################################################CURRENT TEST

#Create dict of motif objects for given list from file
Motifs = {mot:Uniq_Motif(motif=str(mot)) for mot in motif_list}

#Create dict of motifs and their locations for each gene with structure of {Gene1 :  {motif1:[loc1, loc2, loc3], motif2:[loc1,loc2]}}
gene_motif_dict = defaultdict(dict)
for gene,obj in Genes.items():
	for motifobj in Motifs.values(): 
		gene_motif_dict[gene].update(motifobj.gene_search(obj))     #updates current gene key with next motif dict

print(gene_motif_dict)
print('\n')

class Drawing:
	"""
	Drawing class which will take in X and output a 
	pycairo drawing with introns,exons, and motifs labeled. 
	"""
	#takes in gene
	#takes in motif
	#draws 



