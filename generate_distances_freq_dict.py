import sys
import os
import urllib
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def print_results(dict_results, outputname):
    """
    Prints the results with the specified format. Default separation
    is tabs!
    """
    result_ouput = open(outputname, "w")
    for key, value in dict_results.items():
        string = str(key) + "\t"
        for val in value:
            string += str(val) + "\t"
        string += "\n"
        result_ouput.write(string)
    result_ouput.close()

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


def obtain_distances_freq_CIF(filename, outputname):

    transdom_inter_distances = {}
    transdom_aa_freq = {}

    # Select the file and generate a structure var with all the pdb inside

    entity = MMCIFParser()
    structure = entity.get_structure("test", filename)
    model = structure[0]

    transmembrane_dict = determine_transmembrane_domains(filename)

    for key, value in transmembrane_dict.items():
        chain = model[key]

        beg_helix = value[0]
        end_helix = value[1]

        # Iterate each of the resiudes and calculate frequencies
        # and distances between CA atoms

        for i in range(len(beg_helix)):
            hbeg = int(beg_helix[i])
            hend = int(end_helix[i])

            for res1 in range(hbeg, hend + 1):
                res_name1 = chain[res1].get_resname()
                for res2 in range(hbeg, hend + 1):
                    res_name2 = chain[res2].get_resname()

                    if res1 == res2:
                        continue

                    # The class residue automatically calculates the distance
                    dist_value = chain[res1]['CB'] - chain[res2]['CB']

                    # Check if the entry of the 2 aa has been generated reversely!
                    if (res_name2, res_name1) in transdom_inter_distances:
                        transdom_inter_distances[(res_name2, res_name1)].append(dist_value)
                        transdom_aa_freq[res_name1] += 1
                        continue

                    # Save results into a dict-list
                    transdom_inter_distances.setdefault((res_name1, res_name2), [])
                    transdom_inter_distances[(res_name1, res_name2)].append(dist_value)

                    transdom_aa_freq.setdefault(res_name1, 0)
                    transdom_aa_freq[res_name1] += 1



    O_transdom_inter_distances = sorted(transdom_inter_distances.items())
    O_transdom_aa_freq = sorted(transdom_aa_freq.items())

    print_results(transdom_inter_distances, outputname)


if __name__ == "__main__":

    if len(sys.argv) >= 3:
        output = sys.argv[2]
        input_path = sys.argv[1]
        if os.path.isfile(input_path):
            obtain_distances_freq_CIF(input_path, output)
        elif os.path.isdir(input_path):
            dir_to_process = os.getcwd() + '/' + input_path
            for file_dir in os.listdir(dir_to_process):
                obtain_distances_freq_CIF(file_dir, output)
    else:
        raise ValueError("TREX: Please specify output and input")
