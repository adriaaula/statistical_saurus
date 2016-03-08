###we start from a set of hmm output files. Bio.searchIO is a module that works with hmm output files -> no tolera nuestro outputs, asi que 
###sacaremos a la brava

import sys
import os

class Protein(class):
	
	def __init__(self, identifier, profile, evalue): 
		self.__identifier = identifier
		self.__profile = profile
		self.__evalue = evalue
		
	def get_identifier(self): 
		return self.__identifier
		
	def get_profile(self):
		return self.__profile
		
	def get_evalue(self):
		return self.__evalue
		
def create_datasets(): ###crea las carpetas
	
	files = list()
	files = os.listdir("hmmscan")
	for element in files:
		prof = open (element, "r")
		for line in prof:
			if line.startswith("# query sequence file:"):
				l = line
				point = int(l.seqfind("."))
				head = l[point - 5: point - 1]
			if line.startswith("    E-value "):
				next
				next
				l = line
				fields = list()
				fields = l.split() ###split sin argumentos funciona cortando por espacios
				profile = fields[8]
				evalue = fields[0]
		yield Protein(head, profile, evalue)
		os.system("mkdir family/" + profile")
		
for element in create_datasets: ###llena las carpetas
	os.system("cp cif/" + element.get_identifier + ".cif family/" + profile + "/")			
				
				
				
				
				
