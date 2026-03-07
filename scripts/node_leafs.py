#output the count and percentage of strains which carries the ancestor mutatio
#-n on a specific branch

from ete3 import Tree
import sys
import os
import re

lowebr = []
with open('../../data/RLC_lowmapK50E4_H37Rv_pos.txt') as f:
    for line in f:
        lowebr.append(line.strip())

lowebr = set(lowebr)

nwk = open(sys.argv[1]) 
tree = Tree(nwk.read(),format=1)

with open(sys.argv[2]) as f: #db_mutation.txt
    for line in f:
        rows = line.strip().split()
        node = rows[0]
        pos = rows[2]
        pos = pos.split('_')[0]
        if (pos not in lowebr): 
            query = tree.search_nodes(name=node)
            leafs = []
            list(map(leafs.extend, [mnode.get_leaf_names() \
            for mnode in query]))
            leafsnum = len(leafs)
            print('{}\t{}'.format(line.strip(),leafsnum))

