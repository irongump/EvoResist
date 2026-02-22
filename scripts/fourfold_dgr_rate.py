import os
import glob
import pandas as pd
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from collections import Counter

# ==========================================
# 1. Core Parameters and Global Variables
# ==========================================
# Directory containing SNP annotation files
INPUT_DIR = "../ano/"
FILE_EXTENSION = "*.ano" 
SEPARATOR = "\t"

# Column names for headerless files
COL_NAMES = [
    'genome_position', 'ref', 'alt', 'codon_position', 
    'ano', 'codon_change', 'locus_tag', 'gene', 'function', 'group'
]

# Only load these 4 columns to save memory
USE_COLS = ['ref', 'alt', 'codon_change', 'gene']

# Drug resistance genes (use set for fast lookup)
DR_GENES = {"rpoB", "katG", "inhA", "pncA", "gyrA", "gyrB", \
    "embB", "ethA", "Rv0678", "aptE", "pepQ", "rrs", "gid", "rpsL", "eis","tlyA"}

# Fourfold degenerate (FFD) codon prefix combinations
FFD_PREFIXES = {"CT", "GT", "TC", "CC", "AC", "GC", "CG", "GG"}

# ==========================================
# 2. Single File Processing Function (Worker)
# ==========================================
def process_single_file(file_path):
    local_counts = {
        'A': Counter(), 'C': Counter(), 
        'G': Counter(), 'T': Counter()
    }
    
    try:
        # Read file: specify no header (header=None) and assign column names
        df = pd.read_csv(
            file_path, 
            sep=SEPARATOR, 
            header=None, 
            names=COL_NAMES, 
            usecols=USE_COLS, 
            dtype=str
        )
        
        # 1. Exclude drug resistance genes (handle NaN)
        df = df[~df['gene'].isin(DR_GENES)]
        
        # 2. Keep only single base substitutions (SNPs)
        df = df[(df['ref'].str.len() == 1) & (df['alt'].str.len() == 1)]
        
        # 3. Filter out rows with empty codon_change
        df = df.dropna(subset=['codon_change'])
        
        # 4. Strict and efficient FFD filtering logic
        # Format: "CAG-CGG" (length must be 7)
        df = df[df['codon_change'].str.len() == 7]
        
        # Extract first two bases of reference and mutant codons
        # C A G - C G G 
        # 0 1 2 3 4 5 6  (slice: [0:2] = ref prefix, [4:6] = alt prefix)
        ref_prefix = df['codon_change'].str[:2]
        alt_prefix = df['codon_change'].str[4:6]
        
        # FFD condition A: first two bases must be one of 8 special amino acid codons
        cond_ffd = ref_prefix.isin(FFD_PREFIXES)
        
        # FFD condition B: first two bases must be identical before and after mutation
        # (meaning the mutation only occurs at the 3rd position, a synonymous change)
        cond_synonymous = (ref_prefix == alt_prefix)
        
        # Apply filter conditions
        df = df[cond_ffd & cond_synonymous]
        
        # 5. Count the final retained true FFD mutations
        for r, a in zip(df['ref'], df['alt']):
            if r in local_counts and a in local_counts:
                local_counts[r][a] += 1
                
    except Exception as e:
        # Skip corrupted files to keep the main process running
        pass
        
    return local_counts

# ==========================================
# 3. Main Program: Parallel Processing and Result Aggregation
# ==========================================
if __name__ == '__main__':
    print(f"1. Scanning directory {INPUT_DIR} for mutation files...")
    all_files = glob.glob(os.path.join(INPUT_DIR, f"**/{FILE_EXTENSION}"), recursive=True)
    total_files = len(all_files)
    print(f"   Found {total_files} files to process.\n")

    if total_files == 0:
        print("No files found. Please check INPUT_DIR and FILE_EXTENSION.")
        exit()

    global_mutation_counts = {
        'A': {'C': 0, 'G': 0, 'T': 0},
        'C': {'A': 0, 'G': 0, 'T': 0},
        'G': {'A': 0, 'C': 0, 'T': 0},
        'T': {'A': 0, 'C': 0, 'G': 0}
    }

    # Reserve 2 CPU cores for system stability
    num_workers = max(1, cpu_count() - 2)
    print(f"2. Starting parallel processing ({num_workers} CPU cores)...")

    with Pool(processes=num_workers) as pool:
        for local_counts in tqdm(pool.imap_unordered(process_single_file, all_files), total=total_files, desc="Processing Samples"):
            for ref_base, alt_counts in local_counts.items():
                for alt_base, count in alt_counts.items():
                    global_mutation_counts[ref_base][alt_base] += count

    print("\n3. All samples processed! Calculating empirical GTR probability matrix...\n")

    # ==========================================
    # 4. Calculate and Output Final Results
    # ==========================================
    print("# ================= GTR Empirical Probability Matrix ================= #")
    for ref_base in ['A', 'C', 'G', 'T']:
        total_mutations_from_ref = sum(global_mutation_counts[ref_base].values())
        
        if total_mutations_from_ref == 0:
            print(f"# Warning: No FFD mutations recorded for base {ref_base}!")
            continue
            
        for alt_base in ['A', 'C', 'G', 'T']:
            if ref_base != alt_base:
                count = global_mutation_counts[ref_base][alt_base]
                freq = count / total_mutations_from_ref
                print(f"# {ref_base} -> {alt_base}: count = {count:<10} frequency = {freq:.12f}")