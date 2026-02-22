import numpy as np
import pandas as pd
from Bio import SeqIO
import multiprocessing as mp
from tqdm import tqdm
import os

# ---------------------------------------------------------
# 1. Initialization and Global Variables Setup
# (Defined globally so child processes can inherit via Copy-on-Write without memory overhead)
# ---------------------------------------------------------
fasta_file = "../../data/tb.ancestor.fasta"
try:
    record = next(SeqIO.parse(fasta_file, "fasta"))
    reference_genome = np.array(list(str(record.seq).upper()))
except FileNotFoundError:
    print(f"Error: Fasta file {fasta_file} not found. Please check the path.")
    exit()

num_positions = len(reference_genome)  
num_mutations = 2345799    
num_simulations = 1000

bases = np.array(["A", "C", "G", "T"])
ref_numeric = np.searchsorted(bases, reference_genome)

# ---------------------------------------------------------
# 2. Empirical GTR Parameterization
# ---------------------------------------------------------
gtr_probs_raw = np.array([
    [0.393015749751, 0.570833434449, 0.036150815800], # A
    [0.165710425873, 0.309783360596, 0.524506213531], # C
    [0.617719955073, 0.243972295021, 0.138307749906], # G
    [0.055698727469, 0.547079292973, 0.397221979558]  # T
])
gtr_probs = gtr_probs_raw / gtr_probs_raw.sum(axis=1)[:, np.newaxis]

mut_counts_ffd = np.array([123621, 491508, 534200, 137364])
ref_counts = np.array([(reference_genome == b).sum() for b in bases])
raw_rates = mut_counts_ffd / ref_counts
base_mutability = raw_rates / np.mean(raw_rates)
site_base_weights = base_mutability[ref_numeric]

# ---------------------------------------------------------
# 3. Generate Global Lambda Array (Executed Once)
# ---------------------------------------------------------
print("Initializing global mutational landscape...")
lambda_mean = num_mutations / num_positions
gamma_alpha = 0.5  

# Generate the background gamma rate for each site
base_lambda_array = np.random.gamma(shape=gamma_alpha, scale=lambda_mean/gamma_alpha, size=num_positions)
# Scale it using the empirical GTR base mutability
lambda_poisson_array = base_lambda_array * site_base_weights

alt = np.array([
    [1, 2, 3],  
    [0, 2, 3],  
    [0, 1, 3],  
    [0, 1, 2]   
])

# ---------------------------------------------------------
# 4. Define Worker Function for Multiprocessing
# ---------------------------------------------------------
def simulate_single_replicate(sim_idx):
    """
    Simulates mutations for a SINGLE replicate. 
    Keeps memory footprint minimal per CPU core.
    """
    # 1. Draw Poisson mutation counts for this specific replicate
    counts = np.random.poisson(lam=lambda_poisson_array)
    
    # 2. Identify mutated positions and expand them based on counts
    nonzero_idx = np.nonzero(counts)[0]
    repeated_pos = np.repeat(nonzero_idx, counts[nonzero_idx])
    
    # 3. Fetch reference bases and determine substitution
    ref_for_event = ref_numeric[repeated_pos]
    
    rand_vals = np.random.rand(len(repeated_pos))
    cum_probs = np.cumsum(gtr_probs[ref_for_event], axis=1)
    rand_choice = (rand_vals[:, None] > cum_probs).sum(axis=1) 
    new_base_numeric = alt[ref_for_event, rand_choice]
    
    # 4. Group by (position, new_base) to count homoplasies within this tree
    events = np.vstack((repeated_pos, new_base_numeric)).T
    unique_events, homoplasy_counts = np.unique(events, axis=0, return_counts=True)
    
    # 5. Return result as a small DataFrame
    df = pd.DataFrame(unique_events, columns=["position", "new_base"])
    df.insert(0, "sim_index", sim_idx) # Insert at the front
    df["simulated_homoplasy_count"] = homoplasy_counts
    
    return df

# ---------------------------------------------------------
# 5. Main Execution: Parallel Processing
# ---------------------------------------------------------
if __name__ == '__main__':
    # Determine number of CPU cores to use (leave 2 for OS stability)
    num_cores = max(1, mp.cpu_count() - 2)
    print(f"Starting parallel simulation of {num_simulations} replicates using {num_cores} CPU cores...")
    
    results_list = []
    
    # Use Pool to distribute replicates across available CPU cores
    with mp.Pool(processes=num_cores) as pool:
        # imap_unordered is highly efficient and pairs perfectly with tqdm
        for df_res in tqdm(pool.imap_unordered(simulate_single_replicate, range(num_simulations)), total=num_simulations, desc="Simulating"):
            results_list.append(df_res)
            
    print("Simulations complete! Aggregating results...")
    
    # ---------------------------------------------------------
    # 6. Build Final Null Distribution and Save
    # ---------------------------------------------------------
    # Concatenate all parallel results into one comprehensive DataFrame
    mutation_df = pd.concat(results_list, ignore_index=True)
    
    # Calculate empirical null distribution
    dist_df = mutation_df.groupby(['sim_index', 'simulated_homoplasy_count']).size().reset_index(name='num_sites')
    dist_df.to_csv("null_mutation_df_GTR.csv", index=False)
    
    summary_dist = dist_df.groupby('simulated_homoplasy_count')['num_sites'].agg(['mean', 'std']).reset_index()
    summary_dist['lower_95'] = np.maximum(0, summary_dist['mean'] - 1.96 * summary_dist['std'])
    summary_dist['upper_95'] = summary_dist['mean'] + 1.96 * summary_dist['std']
    
    summary_dist.to_csv("expected_null_distribution_GTR_Gamma.csv", index=False)
    mutation_df.to_csv("simulated_mutations_raw_GTR_Gamma.csv", index=False)
    
    print("Statistically sound null models saved successfully.")