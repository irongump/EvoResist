import pandas as pd
from Bio import Phylo, SeqIO
import sys
import os
from os.path import basename
def verify_convergence(target_pos, target_alt, target_ref, target_mut_str, 
                       sequences, mutation_df, pos_mapping, 
                       tree, node_lookup, node_leaves):
    """
    Optimized function to verify convergent mutations.
    """
    # 1. Get nodes where the specific mutation occurred
    nodes_with_mut = mutation_df[mutation_df['mutation'] == target_mut_str]['node_id'].tolist()
    if len(nodes_with_mut) < 2:
        return [] # Return empty list if not convergent

    # 2. Extract fasta index (Ensure it is cast to integer)
    if str(target_pos) not in pos_mapping:
        print(f"Warning: Position {target_pos} not found in pos_mapping.")
        return []
        
    fasta_idx = pos_mapping[str(target_pos)]

    results = []
    
    # 3. Iterate through pairs of convergent branches
    for i in range(len(nodes_with_mut)):
        for j in range(i + 1, len(nodes_with_mut)):
            node_a_id = nodes_with_mut[i]
            node_b_id = nodes_with_mut[j]
            
            # Ensure nodes exist in tree
            if node_a_id not in node_lookup or node_b_id not in node_lookup:
                continue

            node_a = node_lookup[node_a_id]
            node_b = node_lookup[node_b_id]
            
            # Find MRCA
            mrca = tree.common_ancestor(node_a, node_b)
            
            # Identify "Intermediate" samples using pre-cached leaves (O(1) time)
            all_mrca_leaves = node_leaves.get(mrca.name, set())
            leaves_a = node_leaves.get(node_a_id, set())
            leaves_b = node_leaves.get(node_b_id, set())
            
            intermediate_samples = all_mrca_leaves - (leaves_a | leaves_b)
            
            # 4. Check sequences for these intermediate samples
            check_data = {
                "mutation": target_mut_str,
                "pair": f"{node_a_id}-{node_b_id}",
                "mrca": mrca.name if mrca.name else "Unnamed_Node", # Added MRCA node ID
                "WT_base": target_ref,                              # Added Wild Type base
                "total_intermediate": len(intermediate_samples),
                "wt_count": 0,
                "n_count": 0,
                "other_mut_count": 0
            }
            
            for sample in intermediate_samples:
                if sample not in sequences:
                    continue 
                    
                base = str(sequences[sample].seq[fasta_idx]).upper()
                if base == 'N':
                    check_data["n_count"] += 1
                elif base == target_ref:
                    check_data["wt_count"] += 1
                else:
                    check_data["other_mut_count"] += 1
            
            results.append(check_data)
            
    return results

# ==========================================
# Data Loading and Execution Section
# ==========================================

def main():
    if len(sys.argv) < 6:
        print("Usage: python script.py <tree.nwk> <mutation.csv> <sample.fa> <pos.mapping> <target_mutations.csv>")
        sys.exit(1)

    tree_file = sys.argv[1] 
    mutation_csv = sys.argv[2] 
    fasta_file = sys.argv[3] 
    pos_mapping_file = sys.argv[4] 
    target_mutation_file = sys.argv[5] 
    
    output_prefix = basename(tree_file).replace(".tre", "")

    print("Loading tree and pre-computing topologies...")
    tree_obj = Phylo.read(tree_file, "newick")
    
    # Global Cache 1: Node lookup dictionary
    node_lookup = {clade.name: clade for clade in tree_obj.find_clades() if clade.name}
    
    # Global Cache 2: Pre-compute leaves for ALL nodes using post-order traversal
    node_leaves = {}
    for clade in tree_obj.find_clades(order='postorder'):
        if clade.is_terminal():
            node_leaves[clade.name] = {clade.name}
        else:
            leaves = set()
            for child in clade.clades:
                if child.name in node_leaves:
                    leaves.update(node_leaves[child.name])
            if clade.name:
                node_leaves[clade.name] = leaves

    print("Loading fasta...")
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    print("Loading mutation dataframe...")
    mut_df = pd.read_csv(mutation_csv, header=None, sep="\t")
    mut_df.columns = ['node_id', 'branch_length', 'mutation', 'db', 'leaves']
    
    print("Loading position mapping...")
    pos_map_dict = {}
    with open(pos_mapping_file) as f:
        for index, line in enumerate(f):
            clean_line = line.strip()
            if clean_line:
                pos_map_dict[clean_line] = index

    print("Loading target mutations...")
    target_mutation_df = pd.read_csv(target_mutation_file, sep="\t")
    
    print("Verifying convergence...")
    all_results = []
    
    for index, row in target_mutation_df.iterrows():
        pos = str(row['pos'])
        alt = str(row['alt'])
        
        # Robustness: Handle potential NaN values in 'ref' column gracefully and ensure uppercase
        if 'ref' in row and pd.notna(row['ref']):
            ref = str(row['ref']).upper() 
        else:
            ref = ""
            
        target_mut_str = f"{pos}_{alt}"
        
        # Verify and extend results list
        res_list = verify_convergence(
            target_pos=pos,
            target_alt=alt,
            target_ref=ref,
            target_mut_str=target_mut_str,
            sequences=seq_dict, 
            mutation_df=mut_df, 
            pos_mapping=pos_map_dict,
            tree=tree_obj,
            node_lookup=node_lookup,
            node_leaves=node_leaves
        )
        all_results.extend(res_list)
    
    print("Saving results...")
    if all_results:
        final_df = pd.DataFrame(all_results)
        output_file = f"{output_prefix}_convergence_verification_result.csv"
        final_df.to_csv(output_file, index=False)
        print(f"Done! Successfully saved {len(final_df)} comparisons.")
    else:
        print("Done! No convergent mutations met the criteria. No output file generated.")

if __name__ == "__main__":
    main()