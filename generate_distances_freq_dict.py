import sys
import urllib
import logging
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBIO import PDBIO
from numpy import *


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
        chain_helix_domains[chaino[c]].append((beg_helix[c], end_helix[c]))

    logging.info('There are {} chains with transmembrane dom: {}'.format(len(chain_helix_domains),
                                                                         chain_helix_domains.keys()))
    return chain_helix_domains


def calculate_distance(res1, res2, chain):

    res_name1 = chain[res1].get_resname()
    res_name2 = chain[res2].get_resname()

    if res_name1 == 'GLY':
        # get atom coordinates as vectors
        n = chain[res1]['N'].get_vector()
        c = chain[res1]['C'].get_vector()
        ca = chain[res1]['CA'].get_vector()
        # center at origin
        n = n - ca
        c = c - ca
        # find rotation matrix that rotates n -120 degrees along the ca-c vector
        rot = rotaxis(-pi * 120.0 / 180.0, c)
        # apply rotation to ca-n vector
        cb_at_origin = n.left_multiply(rot)
        # put on top of ca atom
        cb1 = cb_at_origin + ca
    if res_name2 == 'GLY':
        # get atom coordinates as vectors
        n = chain[res2]['N'].get_vector()
        c = chain[res2]['C'].get_vector()
        ca = chain[res2]['CA'].get_vector()
        # center at origin
        n = n - ca
        c = c - ca
        # find rotation matrix that rotates n -120 degrees along the ca-c vector
        rot = rotaxis(-pi * 120.0 / 180.0, c)
        # apply rotation to ca-n vector
        cb_at_origin = n.left_multiply(rot)
        # put on top of ca atom
        cb2 = cb_at_origin + ca
        if res_name1 == 'GLY':
            distance = sqrt((cb1[0]-cb2[0])**2 + (cb1[1]-cb2[1])**2 + (cb1[2]-cb2[2])**2)

        else:
            cb1 = chain[res1]['CB'].get_vector()
            distance = sqrt((cb1[0]-cb2[0])**2 + (cb1[1]-cb2[1])**2 + (cb1[2]-cb2[2])**2)
    elif res_name1 == 'GLY':
        cb2 = chain[res2]['CB'].get_vector()
        distance = sqrt((cb1[0]-cb2[0])**2 + (cb1[1]-cb2[1])**2 + (cb1[2]-cb2[2])**2)
    else:
        distance = chain[res1]['CB'] - chain[res2]['CB']

    return distance


def save_results_dict(res_name1, res_name2, distance, dict_name_dist, dict_name_freq):
    """
    Saves the results of the distances and the frequency in the designed dict.
    The frequency dict is given by default, not necessary to definite it
    """
    # Check if the entry of the 2 aa has been generated reversely!
    if (res_name2, res_name1) in dict_name_dist:
        dict_name_dist[(res_name2, res_name1)].append(distance)
        return
    # Save results into a dict-list
    dict_name_dist.setdefault((res_name1, res_name2), [])
    dict_name_dist[(res_name1, res_name2)].append(distance)
    dict_name_freq.setdefault(res_name1, 0)
    dict_name_freq[res_name1] += 1


def print_results(dict_results, outputname):
    """
    Prints the results with the specified format. Default separation
    is tabs!
    """
    result_ouput = open(outputname, "w")

    for key, value in sorted(dict_results.items()):
        string = str(key) + "\t"
        for val in value:
            string += str(val) + "\t"
        string += "\n"
        result_ouput.write(string)
    result_ouput.close()

def obtain_distances_freq_CIF(filename, outputname, outputname2):

    transdom_inter_distances = {}
    transdom_intra_distances = {}
    transdom_aa_freq = {}

    # Select the file and generate a structure var with all the pdb inside
    mmcif_dict = MMCIF2Dict(filename)
    entity = MMCIFParser()
    structure = entity.get_structure(mmcif_dict['_entry.id'], filename)
    logging.info('Working with {}'.format(mmcif_dict['_entry.id']))
    model = structure[0]

    transmembrane_dict = determine_transmembrane_domains(filename)

    for key, value in transmembrane_dict.items():
        logging.info('Chain processed: {}\n Transdom: {}'.format(key, value))
        chain1 = model[key]
        for val in value:

            beg_helix1, end_helix1 = int(val[0]), int(val[1])
            for res1 in range(beg_helix1, end_helix1 + 1):
                for res2 in range(beg_helix1, end_helix1 + 1):
                    if res1 == res2:
                        continue
                    res_name1 = chain1[res1].get_resname()
                    res_name2 = chain1[res2].get_resname()
                    distance = calculate_distance(res1, res2, chain1)
                    save_results_dict(res_name1, res_name2,
                                      distance, transdom_intra_distances,
                                      transdom_aa_freq)
            for val2 in value:
                if val == val2:
                    continue
                beg_helix2, end_helix2 = int(val2[0]), int(val2[1])
                for res1 in range(beg_helix1, end_helix1 + 1):
                    for res2 in range(beg_helix2, end_helix2 + 1):
                        res_name1 = chain1[res1].get_resname()
                        res_name2 = chain1[res2].get_resname()
                        distance = calculate_distance(res1, res2, chain1)
                        save_results_dict(  res_name1, res_name2,
                                            distance, transdom_inter_distances,
                                            transdom_aa_freq)

            print(val)


        # ################INTRA DISTANCES ####################################
        # Iterate each of the residues and calculate frequencies
        # and distances between CA atoms

        #########################################################################
        # #######################INTER DISTANCES  ##############################
        for key, value in transmembrane_dict.items():
            chain2 = model[key]
        # ############################## SAME CHAIN#############################
            if chain1 == chain2:
                continue
            for val in value:
                beg_helix1, end_helix1 = int(val[0]), int(val[1])
                for val in value:
                    pass

    print_results(transdom_intra_distances, outputname)
    print_results(transdom_inter_distances, outputname2)


if __name__ == "__main__":
    import os

    logging.basicConfig(filename=("log_record_" + sys.argv[0]), level=logging.INFO)
    logging.info('\n\nNew running')


    if len(sys.argv) == 4:
        output2 = sys.argv[3]
        output = sys.argv[2]
        input_path = sys.argv[1]
        if os.path.isfile(input_path):
            logging.info('Processing a file')
            obtain_distances_freq_CIF(input_path, output, output2)
        elif os.path.isdir(input_path):
            logging.info('Processing a dir')
            dir_to_process = os.getcwd() + '/' + input_path
            for file_dir in os.listdir(dir_to_process):
                obtain_distances_freq_CIF(file_dir, output, output2)
    else:
        raise ValueError("TREX: Please specify input [1], output (intra_distances [2] and inter_distances [3]) ")
