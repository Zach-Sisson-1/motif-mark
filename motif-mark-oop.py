#!/usr/bin/env python

import argparse
from distutils.command.install_lib import PYTHON_SOURCE_EXTENSION
import re
from urllib.request import parse_http_list

def get_args():
	parser = argparse.ArgumentParser(description="XXXX")
	parser.add_argument("-f", "--fasta_file", help = "Path to the input Fasta File", type=str, required=True)
	parser.add_argument("-m", "--motifs_file", help = "Path to the input Motifs file", type=str, required=True)
	return parser.parse_args()

args = get_args()


#Initalize class for Motifs:

class motif:
	"""
	Motif class that takes in motif pattern and stores all possible combinations as a list
	"""
	def __init__(self, motif_str):
		self.motif = motif_str.lower()
		self.purines = ["a","u","g"] #Purines 
		self.pyrimidines = ["c","t"] #pyrimadines 
		self.motif_list_of_lists = [] # list of lists that will be used to generate all possible motif combinations
		self.motif_combinations_pre = [] 
		self.motif_combinations = [] # list that will contain all possible nucleotide sequences given motif
		
		#Read each nucleotide in the motif and append the appropriate list of possible nucleotides
		for nucleotide in self.motif:
			if nucleotide == "r":
				self.motif_list_of_lists.append(self.purines)
			elif nucleotide =="y":
				self.motif_list_of_lists.append(self.pyrimidines)
			else:
				self.motif_list_of_lists.append(nucleotide)

		#iterate over list of lists, creating a master list of all possible nucleotide sequences, given the motif
		self.motif_combinations_pre = list(itertools.product(*self.motif_list_of_lists))

		#Join together the dictionaries, creating the strings of nucelotides
		for i in self.motif_combinations_pre:
			string = ""
			self.motif_combinations.append(string.join(i))

	def gene_search(self,Gene):
		"""
		Takes in a gene object and searches for motifs in the gene, returning three dictionaries, 
		each containing the unique motifs as keys and the matching locations as values.
		""" 
		#initialize dicts which will house motifs as keys and location hits as values
		pre_exon_dict = {} 
		exon_dict = {}
		post_exon_dict = {}
		
		pre_exon = Gene.pre_exon
		exon = Gene.exon
		post_exon = Gene.post_exon

		#Create iterobjects for each motif combination which locates all the motif hits for the given sequence
		for seq in self.motif_combinations:
			pre_exon_hits = re.finditer(seq,pre_exon) 
			exon_hits = re.finditer(seq,exon)
			post_exon_hits = re.finditer(seq,post_exon)

			#Next add motif with its hits to dict, or "n/a" if none found
			if len(pre_exon_hits) != 0:
				pre_exon_dict[seq] = pre_exon_hits
			else:
				pre_exon_dict[seq] = "n/a"
			if len(exon_hits) != 0:
				exon_dict[seq] = exon_hits
			else:
				exon_dict[seq] = "n/a"
			if len(post_exon_hits) != 0:
				post_exon_dict[seq] = post_exon_hits
			else:
				post_exon_dict[seq] = "n/a"
		return pre_exon_dict, exon_dict, post_exon_dict


class Gene:
	"""
	Gene class that takes in a gene pattern and stores all exons, pre-exons and post-exons.
	"""
	def __init__(self,gene_seq):
		self.gene_seq = gene_seq
		self.pre_exon = ""  
		self.exon = ""
		self.post_exon = ""	

	#Parse gene_seq into pre-exon, exon, and postexon
		self.pre_exon = re.findall('[a-z]+',self.gene_seq)[0]
		self.post_exon = re.findall('[a-z]+',self.gene_seq)[1]
		self.exon = re.findall('[A-Z]+',self.gene_seq)



#read in motif file, creating a dict to hold motif objects
motif_dict = {} #dict to hold input motifs and their class objects
with open(args.motifs_file, "r") as motif_file:
	for line in motif_file:
		Motif = line.strip('\n')
		motif_dict[Motif] = motif(str(Motif))

####working on this#########
#Read in fasta file, creating a dict to hold gene objects
with open(args.fasta_file, "r") as fasta_file:
	gene_dict = {} #empty dict that will contain fasta headers and corresponding genes
	working_seq =""

	for line in fasta_file:
		if line[0] == '>':
			new_key = line[]
			gene_dict[old_key] = working_seq 
			old_key = line[1:]      
			working_seq = ""

		else:
			working_seq += line.strip('\n')

for key,value in gene_dict.items():
	print(key)
	print('\n')
	print(value)
	print('\n')

		

class Drawing:
	"""
	Drawing class which will take in X and output a 
	pycairo drawing with introns,exons, and motifs labeled. 
	"""
	#takes in gene
	#takes in motif
	#draws 



	


#Still need to parse Fafsa file
#Still need to create drawing class that inputs something and creates drawings for pycairo.
	





#Questions:: account for reverse complement in title?
