import re
import sys

# Get ancestor nodes mutations
# Usage: python get_anc_mut.py [candidate_nodes_file] [db_mutation_file]

def get_anc_nodes(candidate_file):
    anc_nodes = []
    with open(candidate_file) as f:
        for line in f:
            if re.search('node_ID',line):
                continue
            nodeid = line.strip().split('\t')[0]
            anc_nodes.append(nodeid)
    return anc_nodes

def get_anc_mut(anc_nodes, db_mutation_file):
    with open(db_mutation_file) as f:
        for line in f:
            node = line.strip().split('\t')[0]
            if node in anc_nodes:
                print(line.strip())
            

if __name__ == '__main__':
    candidate_file = sys.argv[1] if len(sys.argv) > 1 else 'data/100k_5000_level1_candidate_nodes.txt'
    db_mutation_file = sys.argv[2] if len(sys.argv) > 2 else 'data/100k_5000_db_mutation2.txt'
    anc_nodes = get_anc_nodes(candidate_file)
    get_anc_mut(anc_nodes, db_mutation_file)