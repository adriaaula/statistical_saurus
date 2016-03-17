import sys
import urllib
import os
import re
import zipfile
import time
from time import gmtime, strftime
import logging
import pickle
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from numpy import *


def calculate_distance(res1, res2, chain, chain2):
    """
    Retrieves a distance value between two residues.
    Does the calculation for GLY residues.
    """

    res_name1 = chain[res1].get_resname()
    res_name2 = chain2[res2].get_resname()
    try:
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
    except:
        logging.info('Some trouble to access the atoms, located in {} or {}. Most probable explanation: \
        The structure presents only CA'.format(res1,res2))
        return 0

    return distance


def select_and_save_distances(beg_helix1,end_helix1,
                              beg_helix2,end_helix2,
                              chain1,chain2, transdom_intra_distances,
                              transdom_aa_freq, transdom_dist_freq):
    # Compares the positions of each chain and within.
    # Checks if the residue has name, some of them are empty
    if chain1 is chain2:
        if beg_helix1 is beg_helix2:
            beg_helix2 += 5

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
                print('This position ( {} ) didnt present name, skiping it '.format(res2))
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
        dict_name_freq[res_name1]  += 1
        return

    # Save results into a dict-list
    dict_name_dist.setdefault((res_name1, res_name2), {})
    dict_name_dist[(res_name1, res_name2)].setdefault(distance, 0)
    dict_name_dist[(res_name1, res_name2)][distance] += 1
    dict_name_dist[(res_name1, res_name2)].setdefault(999, 0)
    dict_name_dist[(res_name1, res_name2)][999]+= 1

    dict_name_freq.setdefault(res_name1, 0)
    dict_name_freq[res_name1] += 1
    dict_freq_dist.setdefault(distance,0)
    dict_freq_dist[distance] += 1
    return


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

    list_of_dict = list(tuple_dict)

    transdom_intra_distances = list_of_dict[0]
    for dires,value in transdom_intra_distances.items():
        for dist,count in value.items():
            if dist == 999:
                continue
            transdom_intra_distances[dires][dist] = transdom_intra_distances[dires][dist] / transdom_intra_distances[dires][999]

    list_of_dict[0] = transdom_intra_distances
    tuple_dict = list_of_dict


    pickle.dump(tuple_dict, open( "dicts_baby.p", "wb" ) )
    logging.info('THE END\n\n')
    exit()


def MAIN(directory):

    """
    Calculates intra and inter distances of transmembrane domains present in
    the protein.
    Returns dataset of distances for all the proteins analysed.
    """
    transdom_inter_distances = {}
    transdom_intra_distances = {}
    transdom_aa_freq = {}
    transdom_dist_intra_freq = {}
    transdom_dist_inter_freq = {}

    logging.info('Processing the following dir: {}'.format(directory))

    filename_chain_dom_dict = pickle.load( open( "my_phobius_domains.p", "rb" ))
    for PDB_id, chains_file in filename_chain_dom_dict.items():
        if filename_chain_dom_dict[PDB_id] == None:
            continue
        PDB_id = PDB_id[:-5] + 'cif'
        file_dir = os.path.join(directory, PDB_id)
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

        chain_ind = 0
        typical_chains = ('A','B','C','D','E','F','G','H','I','J','K','L')
        for chain, domains in chains_file.items():
            logging.info('Chain processed: {}\n Transdom: {}'.format(chain, domains))
            non_interesting_case = re.search("[0-9]|[a-z]", chain)
            other_non_interesting_case = re.search("[P-Z]", chain)
            if non_interesting_case or other_non_interesting_case:
                continue
            if chain not in model:
                if 'A' not in filename_chain_dom_dict.keys():
                    chain = typical_chains[chain_ind]
                    chain_ind += 1

            chain1 = model[chain]
            for val in domains:
                beg_helix1, end_helix1 = int(val[0]), int(val[1])

                # ################INTRA DISTANCES ####################################

                select_and_save_distances(beg_helix1,end_helix1,
                                          beg_helix1,end_helix1,
                                          chain1,chain1, transdom_intra_distances,
                                          transdom_aa_freq,
                                          transdom_dist_intra_freq)

                # ################INTER DISTANCES SAME CHAIN  ##############################
                for val2 in domains:
                    if val == val2:
                        continue
                    beg_helix2, end_helix2 = int(val2[0]), int(val2[1])
                    select_and_save_distances(beg_helix1,end_helix1,
                                              beg_helix2,end_helix2,
                                              chain1,chain1, transdom_inter_distances,
                                              transdom_aa_freq,
                                              transdom_dist_inter_freq)

                # ################INTER DISTANCES BTWN CHAIN  ########################
                chain_ind2 = 0
                for chain2, value in chains_file.items():
                    non_interesting_case2 = re.search("[0-9]|[a-z]", chain)
                    other_non_interesting_case2 = re.search("[P-Z]", chain)
                    if non_interesting_case2 or other_non_interesting_case2:
                        continue
                    if chain2 not in model:
                        if 'A' not in filename_chain_dom_dict.keys():
                            chain2 = typical_chains[chain_ind]
                            chain_ind2 += 1

                    chain2 = model[chain2]
                    if chain1 == chain2:
                        continue
                    for val3 in value:
                        beg_helix3, end_helix3 = int(val3[0]), int(val3[1])
                        select_and_save_distances(beg_helix1,end_helix1,
                                                  beg_helix3,end_helix3,
                                                  chain1,chain2, transdom_inter_distances,
                                                  transdom_aa_freq,
                                                  transdom_dist_inter_freq)


    tuple_dict = transdom_intra_distances,transdom_inter_distances,transdom_aa_freq,transdom_dist_intra_freq,transdom_dist_inter_freq
    total_division_and_pickling(tuple_dict)


if __name__ == "__main__":


    logging.basicConfig(filename=("log_record_" + sys.argv[0] + ".txt"), level=logging.INFO)
    logging.info('\n\nNew running: ' + strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))
    if len(sys.argv) >= 1:
        dir_path = sys.argv[1]
        if os.path.isdir(dir_path):
            MAIN(dir_path)
        else:
            raise ValueError("We only work with directory files sir.")

    else:
        raise ValueError("TREX: Please specify input [1]")
