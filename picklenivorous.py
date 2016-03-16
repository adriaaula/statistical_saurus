import pickle
from generate_distances_freq_dict import *
from membranosaurus import *
import sys
import urllib
import os
import zipfile
import time
import pexpect
from time import gmtime, strftime
import logging
import pickle
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from numpy import *


def select_and_save_distances_query(beg_helix1,end_helix1,
                              beg_helix2,end_helix2,
                              chain1,chain2, transdom_x_distances):
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

            save_results_query(res1, res_name1, res_name2,
                              distance, transdom_x_distances,
                              transdom_aa_freq,
                              transdom_dist_freq)


def save_results_query(res1, res_name1, res_name2, distance, dict_name_dist):
    """
    Saves the results of the distances and the frequency in the designed dict.
    The frequency dict is given by default, not necessary to definite it
    """
    # Check if the entry of the 2 aa has been generated reversely!
    if (res_name2, res_name1) in dict_name_dist[res1]:
        dict_name_dist[res1][(res_name2, res_name1)].setdefault(distance, 0)
        dict_name_dist[res1][(res_name2, res_name1)][distance] += 1
        dict_name_dist[res1][(res_name2, res_name1)][999] += 1
        return

    # Save results into a dict-list
    dict_name_dist.setdefault(res1, {})
    dict_name_dist[res1].setdefault((res_name1, res_name2), {})
    dict_name_dist[res1][(res_name1, res_name2)].setdefault(distance, 0)
    dict_name_dist[res1][(res_name1, res_name2)][distance] += 1
    dict_name_dist[res1][(res_name1, res_name2)].setdefault(999, 0)
    dict_name_dist[res1][(res_name1, res_name2)][999] += 1
    return



# Recovering of the database with pickle module
tuple_dict = pickle.load( open( "dicts_baby.p", "rb" ))


SP_intra_distances = tuple_dict[0]
SP_inter_distances = tuple_dict[1]
SP_aa_freq = tuple_dict[2]
SP_dist_intra_freq = tuple_dict[3]
#SP_dist_inter_freq = tuple_dict[4]

# Creation of variables to work with new dicts
trans_dom_query = {}
query_inter_distances = {}
query_intra_distances = {}

if len(sys.argv) >= 1:
    dir_path = sys.argv[1]
    if os.path.isdir(dir_path):
        file_dir = os.path.join(directory, file_dir)
        filename = file_dir

    else:
        raise ValueError("We only work with directory files sir.")

else:
    raise ValueError("TREX: Please specify input [1]")



trans_dom_query = determine_transmembrane_domains(filename)

entity = PDBParser()
structure = entity.get_structure('X', filename)
model = structure[0]

for key, value in transmembrane_dict.items():
    logging.info('Chain processed: {}\n Transdom: {}'.format(key, value))
    chain1 = model[key]

    # ################INTRA DISTANCES ####################################

    select_and_save_distances_query(beg_helix1,end_helix1,
                              beg_helix1,end_helix1,
                              chain1,chain1, query_intra_distances)

    # ################INTER DISTANCES SAME CHAIN  ##############################
    for val2 in value:
        if val == val2:
            continue
        beg_helix2, end_helix2 = int(val2[0]), int(val2[1])
        select_and_save_distances_query(beg_helix1,end_helix1,
                                  beg_helix2,end_helix2,
                                  chain1,chain1, query_inter_distances)

    # ################INTER DISTANCES BTWN CHAIN  ########################
    for key, value in transmembrane_dict.items():
        chain2 = model[key]
        if chain1 == chain2:
            continue
        for val3 in value:
            beg_helix3, end_helix3 = int(val3[0]), int(val3[1])
            select_and_save_distances_query(beg_helix1,end_helix1,
                                      beg_helix3,end_helix3,
                                      chain1,chain2, query_inter_distances)


print(query_intra_distances,query_inter_distances)
