import sys
import os
import urllib
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def obtain_distances_freq_CIF(filename, outputname):
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
        uniprot to determine if it is a transmembrane domain or not.
        """
        mmcif_dict = MMCIF2Dict(filename)
        beg_helix = mmcif_dict['_struct_conf.beg_auth_seq_id']
        end_helix = mmcif_dict['_struct_conf.end_auth_seq_id']
        url = 'http://www.uniprot.org/mapping/'

        params = {
        'from':'ACC',
        'to':'P_REFSEQ_AC',
        'format':'tab',
        'query':'P13368 P20806 Q9UM73 P97793 Q17192'
        }

        data = urllib.urlencode(params)
        request = urllib.Request(url, data)
        contact = "" # Please set your email address here to help us debug in case of problems.
        request.add_header('User-Agent', 'Python %s' % contact)
        response = urllib.urlopen(request)
        page = response.read(200000)
        print(page)


    helix_inter_distances = {}
    helix_frequencies_inter_contacts = {}

    # Select the file and generate a structure var with all the pdb inside

    entity = MMCIFParser()
    structure = entity.get_structure("test", filename)
    model = structure[0]
    chain = model['A']

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
                # we erase the HOH residues
                elif res_name1 == 'HOH' or res_name2 == 'HOH':
                    continue

                # The class residue automatically calculates the distance
                dist_value = chain[res1]['CA'] - chain[res2]['CA']

                # Check if the entry of the 2 aa has been generated reversely!
                if (res_name2, res_name1) in helix_inter_distances:
                    helix_inter_distances[(res_name2, res_name1)].append(dist_value)
                    helix_frequencies_inter_contacts[(res_name2, res_name1)] += 1
                    continue

                # Save results into a dict-list
                helix_inter_distances.setdefault((res_name1, res_name2), [])
                helix_inter_distances[(res_name1, res_name2)].append(dist_value)

                helix_frequencies_inter_contacts.setdefault((res_name1, res_name2), 0)
                helix_frequencies_inter_contacts[(res_name1, res_name2)] += 1



    O_helix_inter_distances = sorted(helix_inter_distances.items())
    O_helix_frequencies_inter_contacts = sorted(helix_frequencies_inter_contacts.items())

    print_results(O_helix_inter_distances, outputname)


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
        pass
