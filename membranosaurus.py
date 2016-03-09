###we will start working with a protein whose transmembrane helixes are two:
###
###TM1 : 28 – 50
###
###TM2 : 88 – 111

import sys
import re
import random
import os
import Bio
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

class Membrane(object):
	
	def __init__(self, upper, lower, center): 
		self.__upper = upper
		self.__lower = lower
		self.__center = center
		
	def get_upper(self): 
		return self.__upper
		
	def get_lower(self): 
		return self.__lower
		
	def get_center(self): 
		return self.__center

def get_membrane():
	
	
	parser = MMCIFParser()
	TM = []
	coordinates = []
	vectors = []
	vec = []
	director_vectors = []
	director_vectors_x = []
	director_vectors_y = []
	director_vectors_z = []
	director_vectors_d = []
	files = list()
	files = os.listdir()
	for element in files:
		structure = parser.get_structure('test', element)
		model = structure[0] ### todas las chains son iguales -> no es capaz de diferenciar las cadenas que si que se diferencian en chimera
		mmcif_dict = MMCIF2Dict(element)
		beg_helix = mmcif_dict['_struct_conf.beg_auth_seq_id']
		end_helix = mmcif_dict['_struct_conf.end_auth_seq_id']
		cha_helix = mmcif_dict['_struct_conf.end_auth_asym_id']
		i = 0
		while i < len(beg_helix): ###seleccioinamos las helices transmembrana mediante la longitud (primer filtro)
			if int(end_helix[i]) - int(beg_helix[i]) >= 14:
				a = (beg_helix[i], end_helix[i], cha_helix[i])
				TM.append(a) ###las helices seleccionadas se meten como tuples en una lista
				i += 1
			else:
				i += 1
		#print (TM) ###checkpoint 1 ###las TM seleccionadas son las correctas del modelo
		
		
		for domain in TM:
			
			for res in range(int(domain[0]), int(domain[1])):
				res_coord = list(model[domain[2]][res]['CA'].get_coord())
				coordinates.append(res_coord)
			i = 1
			while i < len(coordinates):
			#print (str(element)) ###Checkpoint 2 -> las listas de coordenadas se crean bien, hay algun error de las funciones própias de biopython
				vec.append(coordinates[i][0] - coordinates[i-1][0])
				vec.append(coordinates[i][1] - coordinates[i-1][1])
				vec.append(coordinates[i][2] - coordinates[i-1][2])
				vectors.append(vec) ###Vectors es la list con los vectores
				i += 1
			for vector in vectors:###Aqui calculo el vector promedio de todos los vectores generados que describen la helix
				director_vectors_x.append(vector[0])
				director_vectors_y.append(vector[1])
				director_vectors_z.append(vector[2])
			director_vectors_d.append(sum(director_vectors_x)/(len(coordinates) -1))
			director_vectors_d.append(sum(director_vectors_y)/(len(coordinates) -1))
			director_vectors_d.append(sum(director_vectors_z)/(len(coordinates) -1))
			director_vectors.append(director_vectors_d)
			coordinates = []
			vec = []
			vectors = []
			director_vectors_d = []
			director_vectors_x = []
			director_vectors_y = []
			director_vectors_z = []
		return director_vectors
			
			
print(get_membrane())
				
				
		
		
