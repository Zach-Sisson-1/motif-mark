#!/usr/bin/env python

import argparse
import itertools

def get_args():
	parser = argparse.ArgumentParser(description="XXXXXXX")
	parser.add_argument("-f", "--fasta_file", help = "Path to the input Fasta File", type=str, required=True)
	parser.add_argument("-m", "--motifs_file", help = "Path to the input Motifs file", type=str, required=True)
	return parser.parse_args()

args = get_args()


#Initalize class for Genes:

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

#read in motif file, creating objects for motif
motif_dict ={} #dict to hold input motifs and their class objects
with open(args.motifs_file, "r") as motif_file:
	for line in motif_file:
		Motif = line.strip('\n')
		motif_dict[Motif] = motif(str(Motif))



#Test
for key,value in motif_dict.items():
	print(value.motif_combinations)


# class Gene:
# 	def __init__(self,intron,exon):
# 		"""Gene class to contain intron and exons """
# 		self.intron = intron
# 		self.	








# #Read in fasta file, making a gene object for each entry
# with open(args.fasta_file, "r") as fasta_file:
# 	gene_dict = {} #empty dict that will contain fasta headers and corresponding genes
# 	for line in fasta_file:
# 		if line[0] == '>':      #new entry = bank the previous sequence as a gene object and reset
# 			gene_dict[]= ""   #set empty string for header value
#  			working_seq = ""
# 		else:
# 			working_seq += line.strip('\n')

		








#Questions:: account for reverse complement in title?
