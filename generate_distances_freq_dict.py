import sys
import urllib
import os
from time import gmtime, strftime
import logging
import pickle
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
    chaino = mmcif_dict['_struct_conf.beg_label_asym_id']

    #Saves each of the transmembrane domains in dict with each chain as key.
    for c in range(len(beg_helix)):
        if int(end_helix[c]) - int(beg_helix[c]) <= 14:
            continue
        chain_helix_domains.setdefault((chaino[c]), [])
        chain_helix_domains[chaino[c]].append((beg_helix[c], end_helix[c]))

    logging.info('There are {} chains with transmembrane dom: {}'.format(len(chain_helix_domains),
                                                                         chain_helix_domains.keys()))
    return chain_helix_domains


def calculate_distance(res1, res2, chain, chain2):
    """
    Retrieves a distance value between two residues.
    Does the calculation for GLY residues.
    """

    res_name1 = chain[res1].get_resname()
    res_name2 = chain2[res2].get_resname()

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
        n = chain2[res2]['N'].get_vector()
        c = chain2[res2]['C'].get_vector()
        ca = chain2[res2]['CA'].get_vector()
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
        cb2 = chain2[res2]['CB'].get_vector()
        distance = sqrt((cb1[0]-cb2[0])**2 + (cb1[1]-cb2[1])**2 + (cb1[2]-cb2[2])**2)
    else:
        distance = chain[res1]['CB'] - chain2[res2]['CB']

    return distance


def select_and_save_distances(beg_helix1,end_helix1,
                              beg_helix2,end_helix2,
                              chain1,chain2, transdom_intra_distances,
                              transdom_aa_freq, transdom_dist_freq):
    # Compares the positions of each chain and within.
    # Checks if the residue has name, some of them are empty
    for res1 in range(beg_helix1, end_helix1 + 1):
        try:
            res_name1 = chain1[res1].get_resname()
        except KeyError as e:
            continue

        for res2 in range(beg_helix2, end_helix2 + 1):
            if res1 == res2:
                continue
            try:
                res_name2 = chain2[res2].get_resname()
            except KeyError as e:
                continue

            distance = calculate_distance(res1, res2, chain1, chain2)
            if distance %1  >= 0.5:
                distance = round(distance)
            else:
                distance = distance - (distance % 1)
                distance = round(distance) + 0.5

            save_results_dict(res_name1, res_name2,
                              distance, transdom_intra_distances,
                              transdom_aa_freq,
                              transdom_dist_freq)


def save_results_dict(res_name1, res_name2, distance, dict_name_dist, dict_name_freq,dict_freq_dist):
    """
    Saves the results of the distances and the frequency in the designed dict.
    The frequency dict is given by default, not necessary to definite it
    """
    # Check if the entry of the 2 aa has been generated reversely!
    if (res_name2, res_name1) in dict_name_dist:
        dict_name_dist[(res_name2, res_name1)].setdefault(distance, 0)
        dict_name_dist[(res_name2, res_name1)][distance] += 1
        dict_name_dist[(res_name2, res_name1)][999] += 1
        dict_name_freq.setdefault(res_name1, 0)
        dict_name_freq[res_name1] += 1
        dict_freq_dist.setdefault(distance,0)
        dict_freq_dist[distance] += 1
        return

    # Save results into a dict-list
    dict_name_dist.setdefault((res_name1, res_name2), {})
    dict_name_dist[(res_name1, res_name2)].setdefault(distance, 0)
    dict_name_dist[(res_name1, res_name2)][distance] += 1
    dict_name_dist[(res_name1, res_name2)].setdefault(999, 0)
    dict_name_dist[(res_name1, res_name2)][999] += 1

    dict_name_freq.setdefault(res_name1, 0)
    dict_name_freq[res_name1] += 1
    dict_freq_dist.setdefault(distance,0)
    dict_freq_dist[distance] += 1
    return

def print_results(dict_results, outputname):
    """
    Prints the results with the specified format. Default separation
    is tabs! Not working right now.
    """
    result_ouput = open(outputname, "w")

    for key, value in sorted(dict_results.items()):
        string = str(key) + "\t"
        for val in value:
            string += str(val) + "\t"
        string += "\n"
        result_ouput.write(string)
    result_ouput.close()

def obtain_distances_freq_CIF(directory):

    """
    Calculates intra and inter distances of transmembrane domains present in
    the protein.
    Returns dataset of distances for all the proteins analysed.
    """
    transdom_inter_distances = {}
    transdom_intra_distances = {}
    transdom_aa_freq = {}
    transdom_dist_freq = {}


    logging.info('Processing the following dir: {}'.format(directory))
    for file_dir in os.listdir(directory):
        file_dir = os.path.join(directory, file_dir)
        filename = file_dir

        # Select the file and generate a structure var with all the pdb inside.
        #SOmetimes there is no option to generate the structure due to PDB files errors.
        mmcif_dict = MMCIF2Dict(filename)
        try:
            entity = MMCIFParser()
            structure = entity.get_structure(mmcif_dict['_entry.id'], filename)
            logging.info('Working with {}'.format(mmcif_dict['_entry.id']))
            model = structure[0]

        except:
            logging.info('NOPE, {}'.format(mmcif_dict['_entry.id']))
            continue

        #CHeck if the protein has helix domains.
        try:
            transmembrane_dict = determine_transmembrane_domains(filename)
        except:
            continue

        for key, value in transmembrane_dict.items():
            logging.info('Chain processed: {}\n Transdom: {}'.format(key, value))
            chain1 = model[key]
            for val in value:
                beg_helix1, end_helix1 = int(val[0]), int(val[1])

                # ################INTRA DISTANCES ####################################
                select_and_save_distances(beg_helix1,end_helix1,
                                          beg_helix1,end_helix1,
                                          chain1,chain1, transdom_intra_distances,
                                          transdom_aa_freq,
                                          transdom_dist_freq)

                # ################INTER DISTANCES SAME CHAIN  ##############################
                for val2 in value:
                    if val == val2:
                        continue
                    beg_helix2, end_helix2 = int(val2[0]), int(val2[1])
                    select_and_save_distances(beg_helix1,end_helix1,
                                              beg_helix2,end_helix2,
                                              chain1,chain1, transdom_inter_distances,
                                              transdom_aa_freq,
                                              transdom_dist_freq)

                # ################INTER DISTANCES BTWN CHAIN  ########################
                for key, value in transmembrane_dict.items():
                    chain2 = model[key]
                    if chain1 == chain2:
                        continue
                    for val3 in value:
                        beg_helix3, end_helix3 = int(val3[0]), int(val3[1])
                        select_and_save_distances(beg_helix1,end_helix1,
                                                  beg_helix3,end_helix3,
                                                  chain1,chain2, transdom_inter_distances,
                                                  transdom_aa_freq,
                                                  transdom_dist_freq)


    tuple_dict = transdom_intra_distances,transdom_inter_distances,transdom_aa_freq,transdom_dist_freq
    total_division_and_pickling(tuple_dict)


def total_division_and_pickling(tuple_dict):
    """Final changes to the dict to have the normalized value division, with probabilities"""

    list_of_dict = list(tuple_dict)

    transdom_inter_distances = list_of_dict[1]
    for dires,value in transdom_inter_distances.items():
        for dist,count in value.items():
            if dist == 999:
                continue
            transdom_inter_distances[dires][dist] = transdom_inter_distances[dires][dist] / transdom_inter_distances[dires][999]

    list_of_dict[1] = transdom_inter_distances
    tuple_dict = list_of_dict

    pickle.dump(tuple_dict, open( "dicts_baby.p", "wb" ) )
    logging.info('THE END\n\n')
    exit()


if __name__ == "__main__":


    logging.basicConfig(filename=("log_record_" + sys.argv[0] + ".txt"), level=logging.INFO)
    logging.info('\n\nNew running: ' + strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))
    if len(sys.argv) >= 1:
        dir_path = sys.argv[1]
        if os.path.isdir(dir_path):
            obtain_distances_freq_CIF(dir_path)
        else:
            raise ValueError("We only work with directory files sir.")

    else:
        raise ValueError("TREX: Please specify input [1]")
