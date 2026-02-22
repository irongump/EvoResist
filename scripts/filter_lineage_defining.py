import sys
import re
from collections import defaultdict

# Filter mutations located in lineage-defining nodes
# based on specific lineages
# Usage: python filter_lineage_defining.py <sublineage_file> <subl_name> \
#   [candidate_nodes_file] [db_mutation_file] [lineage_ancestor_file]

def anc2node(candidate_file):
    anc2node_dict = {}
    with open(candidate_file) as f:
        for line in f:
            if re.search('node_ID',line):
                continue
            node_id = line.split('\t')[0]
            anc = line.split('\t')[1]
            anc2node_dict[anc] = node_id
    return anc2node_dict

def read_anc_mut(db_mutation_file):
    node2mut = defaultdict(list)
    with open(db_mutation_file) as f:
        for line in f:
            node_id = line.strip().split('\t')[0]
            mut = line.strip().split('\t')[2]
            node2mut[node_id].append(mut)
    return node2mut

def read_lineage_anc(lineage_ancestor_file):
    lineage2anc = defaultdict(list)
    with open(lineage_ancestor_file) as f:
        for line in f:
            if re.search('ancestor_nodes',line):
                continue
            cols = line.strip().split(',')
            lineage = cols[0]
            for x in cols[1:]:
                if x and x.strip():
                    lineage2anc[lineage].append(x)
    return lineage2anc

def filter_lineage_defining(sublineage_file, subl_name,
                            candidate_file, db_mutation_file,
                            lineage_ancestor_file):
    anc2node_dict = anc2node(candidate_file)
    node2mut = read_anc_mut(db_mutation_file)
    lineage2anc = read_lineage_anc(lineage_ancestor_file)
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
                else:  # filter mutations in Node1
                    if mut not in ancestor_muts:
                        print(line.strip())

if __name__ == '__main__':
    sublineage_file = sys.argv[1]
    subl_name = sys.argv[2]
    candidate_file = sys.argv[3] if len(sys.argv) > 3 else 'data/100k_5000_level1_candidate_nodes.txt'
    db_mutation_file = sys.argv[4] if len(sys.argv) > 4 else 'data/100k_5000_db_mutation2.txt'
    lineage_ancestor_file = sys.argv[5] if len(sys.argv) > 5 else 'data/lineage_ancestor_nodes.csv'
    filter_lineage_defining(sublineage_file, subl_name,
                            candidate_file, db_mutation_file,
                            lineage_ancestor_file)