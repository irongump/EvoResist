import re
import sys

#get ancestor nodes mutations

def get_anc_nodes():
    anc_nodes = []
    with open('../data/100k_5000_level1_candidate_nodes.txt') as f:
        for line in f:
            if re.search('node_ID',line):
                continue
            nodeid = line.strip().split('\t')[0]
            anc_nodes.append(nodeid)
    return anc_nodes

def get_anc_mut():
    with open('../data/100k_5000_db_mutation2.txt') as f:
        for line in f:
            node = line.strip().split('\t')[0]
            if node in anc_nodes:
                print(line.strip())
            

if __name__ == '__main__':
    anc_nodes = get_anc_nodes()
    get_anc_mut()