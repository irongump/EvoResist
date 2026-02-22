import sys
import re
from collections import defaultdict

#filter mutations which located in lineage defining nodes
#based on specific lineages

def anc2node():
    anc2node_dict = {}
    with open('100k_5000_level1_candidate_nodes.txt') as f:
        for line in f:
            if re.search('node_ID',line):
                continue
            node_id = line.split('\t')[0]
            anc = line.split('\t')[1]
            anc2node_dict[anc] = node_id
    return anc2node_dict

def read_anc_mut():
    node2mut = defaultdict(list)
    with open('100k_5000_db_mutation2.txt') as f:
        for line in f:
            node_id = line.strip().split('\t')[0]
            mut = line.strip().split('\t')[2]
            node2mut[node_id].append(mut)
    return node2mut

def read_lineage_anc():
    lineage2anc = defaultdict(list)
    with open('lineage_ancestor_nodes.csv') as f:
        for line in f:
            if re.search('ancestor_nodes',line):
                continue
            cols = line.strip().split(',')
            lineage = cols[0]
            for x in cols[1:]:
                if x and x.strip():
                    lineage2anc[lineage].append(x)
    return lineage2anc

def filter_lineage_defining(sublineage_file,subl_name):
    anc2node_dict = anc2node()
    node2mut = read_anc_mut()
    lineage2anc = read_lineage_anc()
    ancestors = lineage2anc[subl_name]
    ancestor_nodes = []
    for x in ancestors:
        ancestor_nodes.append(anc2node_dict[x])
    
    ancestor_muts = []
    for y in ancestor_nodes:
        ancestor_muts.extend(node2mut[y])

    with open(sublineage_file) as f:
        for line in f:
            mut = line.strip().split('\t')[2]
            node = line.strip().split('\t')[0]
            if subl_name == 'Lineage5':
                if mut not in ancestor_muts:
                    print(line.strip())
            else:
                if node != 'Node1':
                    print(line.strip())
                else:  #filter mutations in Node1
                    if mut not in ancestor_muts:
                        print(line.strip())

if __name__ == '__main__':
    sublineage_file = sys.argv[1]
    subl_name = sys.argv[2]
    filter_lineage_defining(sublineage_file,subl_name)