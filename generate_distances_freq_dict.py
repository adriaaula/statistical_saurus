import sys
import os
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict




def calculate_distances_frequencies_from_CIF(filename):

    helix_inter_distances = {}
    helix_frequencies_inter_contacts = {}

    #Select the file and generate a structure var with all the pdb inside

    entity = MMCIFParser()
    structure = entity.get_structure("test", filename)

    model = structure[0]
    chain = model['A']

    #Select from mmcif all the transmembrane domains
    #TO DO: Differenciate between helix and transmembrane domains

    mmcif_dict = MMCIF2Dict(filename)

    beg_helix = mmcif_dict['_struct_conf.beg_auth_seq_id']
    end_helix = mmcif_dict['_struct_conf.end_auth_seq_id']

    #Iterate trough everything for take out each of the resiudes and calculate frequencies
    #and distances between CA atoms


    for i in range(len(beg_helix)):
        hbeg = int(beg_helix[i])
        hend = int(end_helix[i])


        for res1 in range(hbeg,hend+1):
            residue_name1 =chain[res1].get_resname()
            for res2 in range(hbeg,hend+1):
                residue_name2 = chain[res2].get_resname()


                if res1 == res2:
                    continue
                elif residue_name1 == 'HOH' or residue_name2 == 'HOH': #we erase the HOH residues
                    continue

                dist_value = chain[res1]['CA'] - chain[res2]['CA'] #The class residue automatically calculates the distance

                #Check if the entry of the 2 aa has been generated reversely!
                if (residue_name2,residue_name1) in helix_inter_distances:
                    helix_inter_distances[(residue_name2,residue_name1)].append(dist_value)
                    helix_frequencies_inter_contacts[(residue_name2,residue_name1)] += 1
                    continue

                helix_inter_distances.setdefault((residue_name1,residue_name2),[])
                helix_inter_distances[(residue_name1,residue_name2)].append(dist_value)

                helix_frequencies_inter_contacts.setdefault((residue_name1,residue_name2),0)
                helix_frequencies_inter_contacts[(residue_name1,residue_name2)] += 1



    O_helix_inter_distances = sorted(helix_inter_distances.items())
    O_helix_frequencies_inter_contacts = sorted(helix_frequencies_inter_contacts.items())

    result_ouput = open("dataset_inter_distances.tab","w")
    for key,value in O_helix_inter_distances:
        string = str(key) + "\t"
        for val in value:
            string += str(val) + "\t"
        string += "\n"
        result_ouput.write(string)
    result_ouput.close()

    #print(O_helix_frequencies_inter_contacts)

if __name__ == "__main__":

    if len(sys.argv) >= 2:
        input_path = sys.argv[1]
        if os.path.isfile(input_path):
            calculate_distances_frequencies_from_CIF(input_path)
        elif os.path.isdir(input_path):
            dir_to_process = os.getcwd() + '/' +input_path
            for file_dir in os.listdir(dir_to_process):
                calculate_distances_frequencies_from_CIF(file_dir)


######Some trying of the possibilities with working cif files

#for res in res_list:
#    pass

#natoms = 0
#nres = 0
#nchain = 0

#for model in structure.get_list():
#    for chain in model.get_list():
#        nchain += 1
#        for residue in chain.get_list():
#            nres += 1
#            for atom in residue.get_list():
#                natoms += 1
#print(natoms, nres, nchain)
