import sys
import re

lowebr = []
with open('/users/g/m/gmin/mtb_ppe/data/RLC_lowmapK50E4_H37Rv_pos.txt') as f:
    for line in f:
        lowebr.append(line.strip())

lowebr = set(lowebr)

with open(sys.argv[1]) as f: # snp file
    for line in f:
        rows = line.strip().split('\t')
        if (rows[0] not in lowebr):
            print(line.strip())
