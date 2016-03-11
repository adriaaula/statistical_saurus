

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
    coor_x = []
    coor_y = []
    coor_z = []
    D = {}
    b = re.compile(".*\.cif")
    files = list()
    files = os.listdir()
    for element in files:
        if b.match(str(element)):
            entity = MMCIFParser()
            structure = entity.get_structure("test", element)
            model = structure[0] ### todas las chains son iguales -> no es capaz de diferenciar las cadenas que si que se diferencian en chimera
            c = 0
            for item, value in determine_transmembrane_domains(element).items(): #iteramos para cada chain
                for domain in value: #iteramos para cada TM de cada chain
                    res = list(range(int(domain[0]), int(domain[1])))
                    for aa in res:
                        res_coord = list(model[item][aa]['CA'].get_coord())
                        coor_x.append(res_coord[0])
                        coor_y.append(res_coord[1])
                        coor_z.append(res_coord[2])
                        coordinates.append(res_coord)
                        c += 1
                    
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
                #return director_vectors ###Checkpoint 3-> genero el numero de vectores correcto, comprovar que la direccion del vector es correcta
                
                for element in director_vectors:
                    director_vectors_x.append(element[0])
                    director_vectors_y.append(element[1])
                    director_vectors_z.append(element[2])
                final_vector.append(sum(director_vectors_x)/(len(director_vectors)))
                final_vector.append(sum(director_vectors_y)/(len(director_vectors)))
                final_vector.append(sum(director_vectors_z)/(len(director_vectors)))
                #return final_vector ###Checkpoint 4-> genero el vector perpendicular, comprovar que la direccion del vector es correcta
                
                ###Ahora tenemos el vector perpendicular al plano, ahora solo nos falta un punto para que el plano pase por el. Tomaremos el punto central de una de las helices
                point = [float(sum(coor_x)/c) , float(sum(coor_y)/c) , float(sum(coor_z)/c)]
                
                #return point ###Checkpoint 5-> genero un punto correctamente

                
                D["central"] = -(final_vector[0]*point[0] + final_vector[1]*point[1] + final_vector[2]*point[2]) 
                
                norm = (final_vector[0]**2 + final_vector[1]**2 + final_vector[2]**2)**0.5
                unitary_vector = []
                for coord in vector:
                    unitary_vector.append(coord/norm)
                
                p_up = [point[0] + 15*unitary_vector[0], point[1] + 15*unitary_vector[1], point[2] + 15*unitary_vector[2]]
                D["up"] = -(vector[0]*p_up[0] + vector[1]*p_up[1] + vector[2]*p_up[2])
                
                p_up_up = [point[0] + 25*unitary_vector[0], point[1] + 25*unitary_vector[1], point[2] + 25*unitary_vector[2]]
                D["up_up"] = -(vector[0]*p_up_up[0] + vector[1]*p_up_up[1] + vector[2]*p_up_up[2])
                
                p_up_out = [point[0] + 45*unitary_vector[0], point[1] + 45*unitary_vector[1], point[2] + 45*unitary_vector[2]]
                D["up_out"] = -(vector[0]*p_up_up[0] + vector[1]*p_up_up[1] + vector[2]*p_up_up[2])

                p_down = [point[0] - 15*unitary_vector[0], point[1] - 15*unitary_vector[1], point[2] - 15*unitary_vector[2]]
                D["down"] = -(vector[0]*p_down[0] + vector[1]*p_down[1] + vector[2]*p_down[2])
                
                p_down_down = [point[0] - 25*unitary_vector[0], point[1] - 25*unitary_vector[1], point[2] - 25*unitary_vector[2]]
                D["down_down"] = -(vector[0]*p_down_down[0] + vector[1]*p_down_down[1] + vector[2]*p_down_down[2])
                
                p_down_out = [point[0] - 45*unitary_vector[0], point[1] - 45*unitary_vector[1], point[2] - 45*unitary_vector[2]]
                D["down_out"] = -(vector[0]*p_up_up[0] + vector[1]*p_up_up[1] + vector[2]*p_up_up[2])
                
                #plane = [vector[0], vector[1], vector[2], D]
                
                #return D ###Checkpoint 6 -> el programa me devuelve correctamente un dictionary con los valores de D de los diferentes planos
                #return item, value ###Checkpoint 6 -> estoy en un entorno en el que tengo acceso a las cadenas (item) y a los TM (value)
                
                ###Conseguimos el primer residuo de una de las cadenas que contienen residuos transmembrana  
                
                first_list = list(str(tuple(structure[0][item])[0]).split())
                first_num = int(first_list[3][7:])  
                #return first_num ###Checkpoint 7 obtengo el numero del primer residuo del chain, siempre y cuando tenga TM
                
                last = len(tuple(structure[0][item]))
                last_list = list(str(tuple(structure[0][item])[last - 1]).split())
                last_num = int(last_list[3][7:])
                #return last_num ###Checkpoint 8 obtengo el numero del ultimo residuo del chain
                
                loops_ter = [(first_num, int(value[0][0]) - 1), (int(value[-1][-1]) + 1, last_num)]
                loops_inter = []
                i = 1
                while i < len(value):
                    tup = (int(value[i-1][1]) + 1, int(value[i][0]) - 1)
                    loops_inter.append(tup)
                    i += 1
                    
                #return loops_inter ###Checkpoint 9 obtengo los residuos correspondientes a los loops, diferencio entre loops terminales y no terminales
                
                for element in loops_ter:
					
                
                
print(get_membrane())
                
                
        
        
