

import sys
import re
import random
import os
import Bio
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import *
#from generate_distances_freq_dict import *

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
		
def determine_transmembrane_domains(filename):
    """
    Compares the helix domains in PDB with the transmembrane domains in
    uniprot to determine if it is a transmembrane domain or not. TO DO
    Returns a dictionary with the transmembrane domains separated by chain
   """
    chain_helix_domains = {}
    mmcif_dict = MMCIF2Dict(filename)
    beg_helix = mmcif_dict['_struct_conf.beg_auth_seq_id']
    end_helix = mmcif_dict['_struct_conf.end_auth_seq_id']
    chaino = mmcif_dict['_struct_conf.end_auth_asym_id']

    for c in range(len(beg_helix)):
        if int(end_helix[c]) - int(beg_helix[c]) <= 14:
            continue
        chain_helix_domains.setdefault((chaino[c]), [])
        chain_helix_domains[chaino[c]].append((beg_helix[c],end_helix[c]))

    return chain_helix_domains

def get_membrane():
		
	TM = []
	coordinates = []
	vectors = []
	vec = []
	director_vectors = []
	director_vectors_x = []
	director_vectors_y = []
	director_vectors_z = []
	director_vectors_d = []
	final_vector = []
	b = re.compile(".*\.cif")
	files = list()
	files = os.listdir()
	for element in files:
		if b.match(str(element)):
			entity = MMCIFParser()
			structure = entity.get_structure("test", element)
			model = structure[0] ### todas las chains son iguales -> no es capaz de diferenciar las cadenas que si que se diferencian en chimera
			determine_transmembrane_domains(element)
#		mmcif_dict = MMCIF2Dict(element)
#		beg_helix = mmcif_dict['_struct_conf.beg_auth_seq_id']
#		end_helix = mmcif_dict['_struct_conf.end_auth_seq_id']
#		cha_helix = mmcif_dict['_struct_conf.end_auth_asym_id']
#		i = 0
#		while i < len(beg_helix): ###seleccioinamos las helices transmembrana mediante la longitud (primer filtro)
#			if int(end_helix[i]) - int(beg_helix[i]) >= 14:
#				a = (beg_helix[i], end_helix[i], cha_helix[i])
#				TM.append(a) ###las helices seleccionadas se meten como tuples en una lista
#				i += 1
#			else:
#				i += 1
#		print (TM) ###checkpoint 1 ###las TM seleccionadas son las correctas del modelo
		return chain_helix_domains
			
			
			c = 0
			for item in chain_helix_domains: #iteramos para cada chain
				for domain in item: #iteramos para cada TM de cada chain
				
				
					for res in range(int(domain[0]), int(domain[1])):
						res_coord = list(model[item][res]['CA'].get_coord())
						coor_x.append(res_coord[0])
						coor_y.append(res_coord[1])
						coor_z.append(res_coord[2])
						coordinates.append(res_coord)
						c += 1
					i = 1
					for element in coordinates:
						coor_x.append(coordinates)
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
				#return director_vectors ###Checkpoint 3-> genero el numero de vectores correcto, comprovar que la direccion del vector es correcta
				
				for element in director_vectors:
					director_vectors_x.append(element[0])
					director_vectors_y.append(element[1])
					director_vectors_z.append(element[2])
				final_vector.append(sum(director_vectors_x)/(len(director_vectors)))
				final_vector.append(sum(director_vectors_y)/(len(director_vectors)))
				final_vector.append(sum(director_vectors_z)/(len(director_vectors)))
				return final_vector ###Checkpoint 4-> genero el vector perpendicular, comprovar que la direccion del vector es correcta
				
				#hacer comando para ver que el vector se me genera en el SENTIDO correcto
				
	###Ahora tenemos el vector perpendicular al plano, ahora solo nos falta un punto para que el plano pase por el. Tomaremos el punto central de una de las helices
				#point = [sum(coor_x)/c, sum(coor_y)/c, sum(coor_z)/c]
				
				#return point
print(get_membrane())
				
				
		
		
