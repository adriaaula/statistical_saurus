import sys
import os
import zipfile
import time
import pexpect
from time import gmtime, strftime
import logging
import pickle
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from numpy import *

def topconn_run(seq):
    fasta = open('topconn_domains/fasta_PDB.fasta', 'w')
    fasta.write('>FASTA_domains\n')
    fasta.write(seq)
    fasta.close()

    child = pexpect.spawn('python topcons2_wsdl.py -m submit -seq topconn_domains/fasta_PDB.fasta')
    child.expect ('jobid = * ', timeout=120)
    id_job = child.readline()
    id_job = str(id_job[0:-2])
    id_job = id_job[2:-1]

    time.sleep(30)

    x = 0
    while x != 1:
        time.sleep(15)
        child2 = pexpect.spawn('python topcons2_wsdl.py -m get -outpath topconn_domains  -jobid ' + id_job)
        result= child2.readlines()
        result = str(result)
        if 'Wait' in result:
            print('00000hhh damm')
            continue
        else:
            print('ualaaaaaa')
            x = 1

    time.sleep(3)
    match = None
    with zipfile.ZipFile('topconn_domains/' + id_job + '.zip') as dom_file:
        with dom_file.open(id_job + '/query.result.txt') as myfile:
            for line in myfile.readlines():
                line = line.decode('utf8')
                if 'TOPCONS predicted topology' in line:
                    match = 'YEY'
                    continue
                if match ==  'YEY':
                    trans_pred = line
                    print(trans_pred)
                    break
    index = 0
    tm = None
    dom = []
    for aa in trans_pred:
        index += 1
        if aa == 'M' and tm == None:
            dom.append(index)
            tm = 'counting'
        if aa != 'M' and tm == 'counting':
            dom.append(index-1)
            final_dom = dom
            dom = []
            tm = None
            yield tuple(final_dom)


def determine_transmembrane_domains(filename):
    """
    Compares the helix domains in PDB with the transmembrane domains in
    uniprot to determine if it is a transmembrane domain or not.
    """

    chain_trans_dom = {}
    true_chain_trans_dom = {}

    mmcif_dict = MMCIF2Dict(filename)
    chains = mmcif_dict['_entity_poly.pdbx_strand_id']
    seqs = mmcif_dict['_entity_poly.pdbx_seq_one_letter_code_can']
    assembly = mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list']

    # Checks if there is present more than one seq for chain.
    # If there is only one (it is a string), it is processed below and the domains are added to the dict.
    if type(seqs) == list:
        i = 0
        for seq in seqs:
            if seq == None or seq == "\n" or seq == "\t" or seq == "" or seq == "?":
                continue
            else:
                for dom in topconn_run(seq):
                    list_chains = chains[i].split(",")
                    for lett in list_chains:
                        if lett == None or lett == "?":
                            continue
                        chain_trans_dom.setdefault(lett, []).append(dom)
                i += 1

    else:
        if seqs == "?":
            return
        for dom in topconn_run(seqs):
            list_chains = chains.split(",")
            for lett in list_chains:
                if lett == None or lett == "?":
                    continue
                chain_trans_dom.setdefault(lett, []).append(dom)

    if type(assembly) == list:
        for chain in assembly[0]:
            if chain in chain_trans_dom.keys():
                true_chain_trans_dom[chain] = chain_trans_dom[chain]
    else:
        true_chain_trans_dom = chain_trans_dom

    print(true_chain_trans_dom)

    logging.info('There are {} chains with transmembrane dom: {}'.format(len(true_chain_trans_dom),
                                                                         true_chain_trans_dom.keys()))
    return true_chain_trans_dom


if __name__ == '__main__':
    trans_domains_topconn = {}

    directory = sys.argv[1]
    logging.basicConfig(filename=("log_record_" + sys.argv[0] + ".txt"), level=logging.INFO)
    logging.info('\n\nNew running: ' + strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))
    logging.info('Processing the following dir: {}'.format(directory))

    index_dir = 0
    for file_dir in os.listdir(directory):
        index_dir += 1
        if index_dir < 15:
            continue
        file_dir = os.path.join(directory, file_dir)
        filename = file_dir
        mmcif_dict = MMCIF2Dict(filename)
        logging.info('Working with {}'.format(mmcif_dict['_entry.id']))
        trans_domains_topconn[mmcif_dict['_entry.id']] = determine_transmembrane_domains(filename)


pickle.dump(trans_domains_topconn, open( "my_topcons.p", "wb" ) )
print(trans_domains_topconn)
logging.info('THE END\n\n')
exit()
