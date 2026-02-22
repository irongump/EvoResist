import sys
import re

# Usage: python remove_low_ebr.py <low_ebr_file> <snp_file>
# Filter out variants in low-mappability regions.
if len(sys.argv) < 3:
    sys.exit(f"Usage: {sys.argv[0]} <low_ebr_file> <snp_file>")

lowebr = set()
with open(sys.argv[1]) as f:
    for line in f:
        lowebr.add(line.strip())

with open(sys.argv[2]) as f:
    for line in f:
        rows = line.strip().split('\t')
        if rows[0] not in lowebr:
            print(line.strip())
