import pickle
from generate_distances_freq_dict import *
from phobius_tester import *
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
import numpy as np


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
                              distance, transdom_x_distances)


def save_results_query(res1, res_name1, res_name2, distance, dict_name_dist):
    """
    Saves the results of the distances and the frequency in the designed dict.
    The frequency dict is given by default, not necessary to definite it
    """
    # Check if the entry of the 2 aa has been generated reversely!
    dict_name_dist.setdefault(res1, {})
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

def distance_calculator(directory):

    # Creation of variables to work with new dicts
    trans_dom_query = {}
    query_inter_distances = {}
    query_intra_distances = {}

    file_dir = os.path.join(directory, '4mbs.cif')
    filename = file_dir

    trans_dom_query = phobius_runner('/home/adria/UPF/PDB_fasta/4mbs.' + 'fasta')

    entity = MMCIFParser()
    structure = entity.get_structure('X', filename)
    model = structure[0]

    for key, value in trans_dom_query.items():
        logging.info('Chain processed: {}\n Transdom: {}'.format(key, value))
        chain1 = model[key]
        for val in value:
            beg_helix1, end_helix1 = int(val[0]), int(val[1])
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
            for key, value in trans_dom_query.items():
                chain2 = model[key]
                if chain1 == chain2:
                    continue
                for val3 in value:
                    beg_helix3, end_helix3 = int(val3[0]), int(val3[1])
                    select_and_save_distances_query(beg_helix1,end_helix1,
                                              beg_helix3,end_helix3,
                                              chain1,chain2, query_inter_distances)


    return (query_intra_distances,query_inter_distances)

def PMF_pair_calc(res_pair, distance, query_dict, SP_distances_2res, SP_freq_aa, SP_freq_dist):

    PMF_pair = 0
    for dist_query, ocurrences in distance.items():
        acumulative_distance_2res = 0
        acumulative_distance_gen = 0
        if res_pair not in SP_distances_2res.keys():
            res_pair = (res_pair[1],res_pair[0])

        acumulative_distance_2res += sum( freq for dist,freq in  SP_distances_2res[res_pair].items()
                                          if dist <= dist_query and dist != 999 )
        acumulative_distance_gen += sum( freq for dist,freq in  SP_freq_dist.items()
                                          if dist <= dist_query)

        #print(SP_freq_aa[res_pair[0]], SP_freq_aa[res_pair[1]], acumulative_distance_gen)
        PMF_pair_dist = - log(acumulative_distance_2res/
                             (SP_freq_aa[res_pair[0]]*SP_freq_aa[res_pair[1]]*acumulative_distance_gen)
                             )*ocurrences

        PMF_pair += PMF_pair_dist


    return PMF_pair

    #acumulative_distance_2res = np.sum([ SP_distances_2res[val] for val in SP_distances_2res.values() if val <= max_distance ] )

def MAIN(directory):

    # Recovering of the database with pickle module
    tuple_dict = pickle.load( open( "dicts_baby.p", "rb" ))


    SP_intra_distances = tuple_dict[0]
    SP_inter_distances = tuple_dict[1]
    SP_aa_freq = tuple_dict[2]
    total_aa = sum(values for values in SP_aa_freq.values())
    for key,value in SP_aa_freq.items():
        SP_aa_freq[key] =  SP_aa_freq[key] / total_aa

    SP_dist_intra_freq = tuple_dict[3]
    total_intra = sum(values for values in SP_dist_intra_freq.values())
    for key,value in SP_dist_intra_freq.items():
        SP_dist_intra_freq[key] =  SP_dist_intra_freq[key] / total_intra
    SP_dist_inter_freq = tuple_dict[4]
    total_inter = sum(values for values in SP_dist_inter_freq.values())
    for key,value in SP_dist_inter_freq.items():
        SP_dist_inter_freq[key] =  SP_dist_inter_freq[key] / total_inter

    distances_query = distance_calculator(directory)
    query_distances_intra, query_distances_inter = distances_query[0] , distances_query[1]

    PMF_positions = {}

    for pos in query_distances_intra.keys():
        for aa_pairs,distances in query_distances_intra[pos].items():
            PMF_positions[pos] = np.sum(PMF_pair_calc(aa_pairs,distances,
                                                query_distances_intra,
                                                SP_intra_distances,
                                                SP_aa_freq ,
                                                SP_dist_intra_freq) for aa_pairs,distance in query_distances_intra[pos].items())

            PMF_positions[pos] += np.sum(PMF_pair_calc(aa_pairs,distances,
                                                query_distances_inter,
                                                SP_inter_distances,
                                                SP_aa_freq ,
                                                SP_dist_inter_freq) for aa_pairs,distance in query_distances_inter[pos].items())



        print(pos, PMF_positions[pos])



if __name__ == '__main__':


    if len(sys.argv) >= 1:
        dir_path = sys.argv[1]
        if os.path.isdir(dir_path):
            MAIN(dir_path)
        else:
            raise ValueError("We only work with directory files sir.")

    else:
        raise ValueError("TREX: Please specify input [1]")
