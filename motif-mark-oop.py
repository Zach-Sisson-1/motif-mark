#!/usr/bin/env python

import argparse
from collections import defaultdict
import re
import itertools
import cairo
import numpy as np

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
		#update the span to reflect motif length
		index_dict2 ={}
		for motif,spanlist in index_dict.items():
			index_dict2[motif] = []
			for hit in spanlist:
				hit = list(hit)
				hit[1] = int(hit[0] + len(motif))
				hit = tuple(hit)
				index_dict2[motif].append(hit)
		
		return index_dict2
	

class Gene:
	"""
	Gene class that takes in a gene pattern and stores all exons, pre-exons and post-exons.
	"""
	def __init__(self,gene_seq, header):
		self.gene_seq = gene_seq
		self.gene_seq_lower = gene_seq.lower()   #hold lowercase string for motif searching
		self.header = header						#header for FASTA seq 
		self.length = len(self.gene_seq_lower)   

	#Parse gene_seq into pre-exon, exon, and postexon
		self.pre_exon = re.findall('[a-z]+',self.gene_seq)[0]
		self.post_exon = re.findall('[a-z]+',self.gene_seq)[1]
		self.exon = re.findall('[A-Z]+',self.gene_seq)
		

	#calculate lengths which will be used for drawing
		self.pre_exon_length = len(self.pre_exon)
		self.exon_length = len(self.exon[0])
		self.post_exon_length = len(self.post_exon)


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


#Create dict of Gene object for each gene in dictionary.   Keys = header name and values = gene objects 
Genes = {name:Gene(gene_seq=(gene_dict[name]), header=name) for name in gene_dict.keys()}

	
###########################################################CURRENT TEST

#Create dict of motif objects for given list from file
Motifs = {mot:Uniq_Motif(motif=str(mot)) for mot in motif_list}

#Create dict of motifs and their locations for each gene with structure of {Gene1 :  {motif1:[loc1, loc2, loc3], motif2:[loc1,loc2]}}
gene_motif_dict = defaultdict(dict)
for gene,obj in Genes.items():
	for motifobj in Motifs.values(): 
		gene_motif_dict[gene].update(motifobj.gene_search(obj))     #updates current gene key with next motif dict


#Join dictionaries so key = Gene name, value = [Gene object, motif dict]
master_dict = {}
for gene,object in Genes.items():
	master_dict[gene] = [object,gene_motif_dict[gene]]


class Drawing:
	"""
	Drawing class which will take in a dictionary of Genes, their objects, and their motif locations and output a 
	pycairo drawing with introns,exons,and motifs labeled. 
	"""
	def __init__(self,master_dict) -> None:
		#establish the context/size of layout based on number of genes and length of longest gene
		self.master_dict = master_dict				#	store the dict
		counter = 0									#	initialize counter to track longest gene
		for gene,lst in master_dict.items():
			if lst[0].length > counter:
				counter = lst[0].length
		self.X = int(counter+(0.2*counter))								# longest gene + 30%
		self.Y = int(len(master_dict.keys())*150)						# scale height according to # of genes. 
		
		self.spacepergene = (self.Y-20)/(len(master_dict.keys()))        # allot space per gene 
		self.namespace = (self.spacepergene*0.4)						# space between gene name and sequence
		self.betweenspace = (self.spacepergene*0.7)						# space between genes									
		self.exonspace = (self.spacepergene*0.5) 						# height of exon rectangle
		self.motiflen = (self.spacepergene*0.3)							# height of motif rectangle
		self.gapy = 20													# down shift 
		self.gapx = (self.X*0.02) 										# right shift

		#Establish color dict for each motif used for drawing
		np.random.seed(63)
		self.color_dict = {}
		for gene,dict in self.master_dict.values():
			for motif in dict.keys():
				self.color_dict[motif] = tuple(np.random.random(size=3))


	
	def draw(self):	
		#set width and height for image size
		WIDTH, HEIGHT = self.X,self.Y

		#Creates image surface with above dimensions and RGB24 format where each pixel is 32-bit quantity (see https://pycairo.readthedocs.io/en/latest/reference/enums.html#cairo.Format for more info)
		surface = cairo.ImageSurface(cairo.FORMAT_RGB24, WIDTH, HEIGHT)

		#Set drawing context to be the image surface above
		ctx = cairo.Context(surface)
		
		#Draw components for each gene
		yslider = self.gapy							#initialize first y-coordinate
		for gene,lst in master_dict.items():
																#Draw Gene name
			ctx.move_to(self.gapx,yslider)								
			ctx.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
			ctx.set_font_size(15)
			ctx.set_source_rgb(1,1,1) #white color
			ctx.show_text(gene)
			yslider += self.namespace	 #move down for start of gene line
			
																#Draw pre-exon and post-exon Gene line
			#draw pre-exon
			ctx.set_source_rgb((135/255),(131/255),(131/255)) #Gray color
			ctx.set_line_width(5) #set line thickness
			ctx.move_to(self.gapx,yslider)
			ctx.line_to(lst[0].pre_exon_length,yslider) 
			ctx.stroke()
			ctx.fill()
			
			#set rectangle coords for exon	
			hgt = self.exonspace
			wth = lst[0].exon_length
			x0 = lst[0].pre_exon_length
			y0 = yslider-(hgt/2)
			#draw exon rectangle
			ctx.rectangle(x0,y0, wth, hgt)  # Rectangle(x0, y0, width, height)
			ctx.fill()
			
			#draw post-exon						
			ctx.move_to((lst[0].pre_exon_length+lst[0].exon_length),yslider)
			ctx.line_to((lst[0].pre_exon_length+lst[0].exon_length+lst[0].post_exon_length),yslider) 
			ctx.stroke()
			ctx.fill()


		### Draw Motifs for each gene
			for motif,hitlist in lst[1].items():
				ctx.set_source_rgb( (self.color_dict[motif][0]), (self.color_dict[motif][1]), (self.color_dict[motif][2]) ) # Set motif color
				for hitpair in hitlist:
					ctx.rectangle( (hitpair[0]+self.gapx), (yslider-(self.motiflen/2)), (hitpair[1]-hitpair[0]), self.motiflen)
					ctx.fill()

		# move slider to next gene
			yslider += self.betweenspace	
			
	### Draw Legend box
		legend_slider = self.gapy
		ctx.move_to((WIDTH-150),legend_slider)
		ctx.set_source_rgb(1,1,1) #white color
		ctx.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
		ctx.set_font_size(20)
		ctx.show_text("Legend")
		
		legend_slider +=20
		
		ctx.set_source_rgb((135/255),(131/255),(131/255)) #Gray color
		ctx.set_line_width(5) #set line thickness
		ctx.move_to((WIDTH-150),legend_slider)
		ctx.line_to(WIDTH-125,legend_slider) 
		ctx.stroke()
		ctx.fill()
		
		legend_slider +=5

		ctx.move_to((WIDTH-120),legend_slider)
		ctx.show_text("Intron")

		legend_slider += 10
		ctx.rectangle(WIDTH-150,legend_slider, 20, 20)  # Rectangle(x0, y0, width, height)
		ctx.fill()
		
		legend_slider += 15
		
		ctx.move_to((WIDTH-120),legend_slider)
		ctx.set_font_size(20)
		ctx.show_text("Exon")

		for motif,color in self.color_dict.items():
			legend_slider += 15
			ctx.set_source_rgb(color[0],color[1],color[2])
			ctx.rectangle(WIDTH-150,legend_slider, 20, 20)  # Rectangle(x0, y0, width, height)
			ctx.fill()
		
			legend_slider += 15
		
			ctx.move_to((WIDTH-120),legend_slider)
			ctx.set_font_size(20)
			ctx.show_text(motif)
	
		#Save image
		return surface


#Create the object
drawing = Drawing(master_dict=master_dict)

#call the objects method
surface = drawing.draw()

#Specify output file name
filename_output = "{}.png".format(re.match("^([^.]+)",args.fasta_file).group(0))

#save image object
surface.write_to_png(filename_output)

