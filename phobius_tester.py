import sys
import os
import time
import subprocess
from time import gmtime, strftime
import logging
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pickle

def FASTA_iterator(filename):
    """Returns fasta sequence on at a time using yield"""
    fasta_file=open(filename, "r")
    id_fasta=""
    seq_fasta=""
    for line in fasta_file:
        if line.startswith(">"):
            if id_fasta == "":
                id_fasta=line.strip()
                continue

            fasta = id_fasta , seq_fasta
            yield fasta
            seq_fasta=""
            id_fasta=line.strip()

        else:
            seq_fasta += line.strip()

def phobius_runner(filename):
    chain_dom_dict = {}

    for fasta in FASTA_iterator(filename):

        good_entry = 'False'
        fasta_id = fasta[0]
        fasta_chain = fasta_id.split("|")[0].split(":")[-1]
        fasta = fasta[0] + "\n" + fasta[1]

        fasta_file = open('phobius/fasta_PDB.fasta', 'w')
        fasta_file.write(fasta)
        fasta_file.close()

        child = subprocess.check_output("cat phobius/fasta_PDB.fasta | perl phobius/phobius.pl ", shell=True)
        child = child.decode('utf8')
        child = child.split('\n')
        for entry in child:
            if 'TRANSMEM' in entry:
                entry = [en for en in entry.split(" ") if en != ""]
                trans_dom = (entry[2],entry[3])
                chain_dom_dict.setdefault(fasta_chain,[])
                chain_dom_dict[fasta_chain].append(trans_dom)
                good_entry = 'True'

        if good_entry == 'True':
            return chain_dom_dict

if __name__ == '__main__':

    trans_domains_phobius = {}
    directory = sys.argv[1]
    logging.basicConfig(filename=("log_record_" + sys.argv[0] + ".txt"), level=logging.INFO)
    logging.info('\n\nNew running: ' + strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))
    logging.info('Processing the following dir: {}'.format(directory))

    index_dir = 0
    for file_dir in os.listdir(directory):
        index_dir += 1
        if index_dir < 0:
            continue
        file_dir = os.path.join(directory, file_dir)
        filename = file_dir

        mmcif_dict = MMCIF2Dict(filename)
        logging.info('Working with {}'.format(filename))
        trans_domains_phobius[filename.split("/")[-1]] = phobius_runner(filename)

    pickle.dump(trans_domains_phobius, open( "my_phobius_domains.p", "wb" ))
    print('The end motherfuckers')
    exit()
