# General modules
import pickle
import sys
import os
import time
from time import gmtime, strftime
import logging
import argparse
# Graph modules
import plotly
from plotly.graph_objs import Scatter, Layout
# Bio.PDB
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from numpy import *
# Own modules, scripts
from phobius_tester import *
from SP_calculator_DB import calculate_distance


def select_and_save_distances_query(beg_helix1,end_helix1,
                              beg_helix2,end_helix2,
                              chain1,chain2, transdom_x_distances):
    """
        Iterates through all the residues from two domains and saves the distance
        calculated in the calculate_distance function.
    """
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
    Saves the results of the distances of two residues in a given position.
    """
    # Check if the entry of the 2 aa has been generated reversely!
    dict_name_dist.setdefault(res1, {})
    if (res_name2, res_name1) in dict_name_dist[res1]:
        dict_name_dist[res1][(res_name2, res_name1)].setdefault(distance, 0)
        dict_name_dist[res1][(res_name2, res_name1)][distance] += 1
        return

    # Save results into a dict-list
    dict_name_dist.setdefault(res1, {})
    dict_name_dist[res1].setdefault((res_name1, res_name2), {})
    dict_name_dist[res1][(res_name1, res_name2)].setdefault(distance, 0)
    dict_name_dist[res1][(res_name1, res_name2)][distance] += 1
    return

def get_seq_PDB(filename):

    """
    From a PDB created with Modeller or other, iterates trough all atoms and generates
    a fasta sequence with all the information.
    """

    pdb = open (filename, "r")
    seq = ">" + filename[0:4] + ":A|something|other|another\n"
    chain = "A"
    correspondence = { "ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E",
     "PHE":"F", "GLY":"G", "HIS":"H", "ILE":"I", "LYS":"K", "LEU":"L",
     "MET":"M", "ASN":"N", "PRO":"P", "GLN":"Q", "ARG":"R", "SER":"S",
     "THR":"T", "VAL":"V", "TRP":"W", "TYR":"Y"}
    for line in pdb:
        if line.startswith("ATOM"):
            field_list = line.split()

            if field_list[2] == "CA":
                seq += str(correspondence[field_list[3]])
    seq += '\n>\n'
    return seq

def distance_query_generator(filename):
    """
    Creates the following dicts with the values of distances:
        - query_intra_distances: values from the distances same chain same dom.
        - query_inter_distances: values between different chains diff domains.

    """
    # Creation of variables to work with new dicts
    trans_dom_query = {}
    query_inter_distances = {}
    query_intra_distances = {}

    file_id = (filename.split('/')[-1])[:-3]


    if '.pdb' in filename:
        fasta_fil = open('fasta_to_generate.fasta', 'w')
        seq_fasta = get_seq_PDB(filename)
        fasta_fil.write(seq_fasta)
        fasta_fil.close()

        trans_dom_query = phobius_runner('fasta_to_generate.fasta')
        entity = PDBParser()

    elif '.cif' in filename:
        trans_dom_query = phobius_runner('/home/adria/UPF/PDB_fasta/' + file_id + 'fasta')
        entity = MMCIFParser()

    structure = entity.get_structure('X', filename)
    model = structure[0]

    for key, value in trans_dom_query.items():
        logging.info('Chain processed: {}\n Transdom: {}'.format(key, value))
        if '.pdb' in filename:
            for chain in model:
                chain1 = chain
        else:
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
                if '.pdb' in filename:
                    for chain in model:
                        chain2 = chain
                else:
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
    """
    Calculates all the PMF_pair from a given pair of residues.
    In order to do that, calculates:
    - acumulative_distance_2res: frequencies presented at the specific distance or less
      between two residues.
    - acumulative_distance_gen: general frequency of the given distance or less in all the datasets.

    """

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

        PMF_pair_dist = - log(acumulative_distance_2res/
                             (SP_freq_aa[res_pair[0]]*SP_freq_aa[res_pair[1]]*acumulative_distance_gen)
                             )*ocurrences

        PMF_pair += PMF_pair_dist

    return PMF_pair


def apply_window_size(PMF,window):
    """
    Apply an specific windows size selected by the user for showing in the graph the results.
    Calculates the median of the window and assigns at each position of interest the value.
    """
    PMF_median = {}
    list_keys = []
    median_window = 0
    for i in range(0,len(PMF)):
        list_keys = sorted(PMF.keys())
        if len(PMF) >= i+window:
            PMF_median[list_keys[i]] = sum(PMF[list_keys[x]] for x in range(i,i+window)) / window
        else:
            too_much = i+window
            while len(PMF) < too_much:
                too_much -= 1
            PMF_median[list_keys[i]] = sum(PMF[list_keys[x]] for x in range(i,too_much)) / (too_much-i)

    return PMF_median


def make_plot(PMF_positions,file_id,output_path):
    """
    Creates using plotly an offline graph in html with the displaying of the results.
    """
    # Assigns None value to the positions no transmembrane for showing efficiently the results
    for i in range(min(PMF_positions.keys()),max(PMF_positions.keys())):
        if i in PMF_positions.keys():
            continue
        else:
            PMF_positions[i] = None
    positions = [key for key in sorted(PMF_positions.keys())]

    pmf_values = [value for key,value in sorted(PMF_positions.items())]
    color_list = [  'rgb(22, 96, 167)', 'rgb(205, 12, 24)',
                    'rgb(0,128,128)', 'rgb(250,128,114)', 'rgb(154,205,50)'
                    'rgb(255,105,180)', 'rgb(0,191,255)', 'rgb(50,205,50)', 'rgb(165,42,42)']

    # Create a trace

    plotly.offline.plot({
    "data": [
        Scatter(
                x=positions,
                y=pmf_values,
                line = dict(color = random.randint(0,len(color_list)))
                )

    ],
    "layout": Layout(
        xaxis =dict(title = 'Position' ) ,
        yaxis = dict(title = 'Potential Mean Force (kT)'),
        title= '<b>Transmembrane statistic potentials:  ' + file_id + '</b>'
    )
    }, filename = output_path + '.html')


def MAIN(filename,window,output_path):
    """
    Main program with all the subprocesses. Calculates PMF intra and inter
    positions. Creates the plot.
    """
    # Recovering of the database with pickle module
    tuple_dict = pickle.load( open( "SP_dataset.p", "rb" ))

    # Prepares all the dataset for calculate PMF
    SP_intra_distances = tuple_dict[0]
    SP_inter_distances = tuple_dict[1]
    SP_aa_freq = tuple_dict[2]
    SP_dist_intra_freq = tuple_dict[3]
    SP_dist_inter_freq = tuple_dict[4]

    sys.stderr.write('\tCreating the distance query dict\n')
    # Obtain all the distance values
    distances_query = distance_query_generator(filename)
    query_distances_intra, query_distances_inter = distances_query[0] , distances_query[1]

    PMF_positions = {}
    sys.stderr.write('\tCalculating the PMF for each position\n')
    # Calculates all the PMF values for a given position
    for pos in query_distances_intra.keys():
        PMF_positions[pos] = sum(PMF_pair_calc(aa_pairs,distances,
                                               query_distances_intra,
                                               SP_intra_distances,
                                               SP_aa_freq ,
                                               SP_dist_intra_freq) for aa_pairs,distances in query_distances_intra[pos].items())

    for pos in query_distances_inter.keys():
        PMF_positions[pos] += sum(PMF_pair_calc(aa_pairs,distances,
                                            query_distances_inter,
                                            SP_inter_distances,
                                            SP_aa_freq ,
                                            SP_dist_inter_freq) for aa_pairs,distances in query_distances_inter[pos].items())

    sys.stderr.write('\tPlotting the results...')
    #Creates the graph in plotly
    PMF_median= apply_window_size(PMF_positions,window)
    file_id = filename.split('/')[-1]
    make_plot(PMF_median,file_id,output_path)

def argparse_creator():

	parser = argparse.ArgumentParser(description="""This program analizes the folding of transmembrane
                                                    domains in proten models according to knowledge-based potentials""")

	parser.add_argument('-i','--input',
						dest = "infile",
						action = "store",
						default = None,
						help = "Input file or directory. Files must have only one chain")

	parser.add_argument('-o','--output',
						dest = "outfile",
						action = "store",
						default = "output.html",
						help = "Output file or directory")

	parser.add_argument('-w','--window',
						dest = "window",
						action = "store",
						default = 10,
						help = "Size of the window in the plot for statistic potentials")

	options = parser.parse_args()

	infile = options.infile
	output = options.outfile
	window = options.window

	arg = (infile, output, window)
	return arg


if __name__ == '__main__':

    argv = argparse_creator()

    if len(sys.argv) >= 1:
        path = argv[0]
        output_path = argv[1]
        window = int(argv[2])
        if os.path.isdir(path):
            for filename in os.listdir(path):
                sys.stderr.write('Working with {}\n'.format(filename))
                file_pdb = os.path.join(path, filename)
                MAIN(file_pdb,window,output_path)
        elif os.path.isfile(path):
            if '.pdb' in path or '.cif' in path:
                sys.stderr.write('Working with {}\n'.format(path))
                MAIN(path,window,output_path)
            else:
                raise ValueError("This type of file is not accepted")

        else:
            raise ValueError("We only work with directory files sir.")

    else:
        raise ValueError("TREX: Please specify input [1]")
