# ==============================================================================
# Script: Calculate Empirical P-values and FDR for Convergent Mutations
# ==============================================================================
# Install data.table if not already installed: install.packages("data.table")
library(data.table)
library(stringr)
library(dplyr)
# 1. Define file paths (Modify these paths according to your environment)
obs_file <- "../lineage_ann/all_ann_convergent_flt.txt" # Columns expected: pos, ref, alt, count
sim_file <- "simulated_mutations_raw_GTR_Gamma.csv" # Columns: sim_index, position, new_base, simulated_homoplasy_count

# Total number of simulation replicates run in the Python script
num_simulations <- 1000

# ------------------------------------------------------------------------------
# 2. Load data efficiently using fread
# ------------------------------------------------------------------------------
cat("Loading observed and simulated datasets...\n")
obs_dt <- read.delim(obs_file, sep=" ", header=T) %>% 
  filter(!grepl("-", pos) & alt %in% c("A", "C", "G", "T") & ref!=alt) %>%  
  rename(position = pos) %>% mutate(position=as.numeric(position))
obs_dt <- as.data.table(obs_dt)
sim_dt <- fread(sim_file)


# ------------------------------------------------------------------------------
# 3. Align base encoding between Observed and Simulated data
# ------------------------------------------------------------------------------
# In the Python simulation, bases were encoded as integers: A=0, C=1, G=2, T=3
# We need to map the 'alt' character in observed data to these integers
base_map <- c("A" = 0, "C" = 1, "G" = 2, "T" = 3)
obs_dt[, new_base := base_map[alt]]

# ------------------------------------------------------------------------------
# 4. Calculate Empirical P-values (Vectorized & Highly Optimized)
# ------------------------------------------------------------------------------
cat("Calculating site-wise empirical P-values...\n")

# Set keys for ultra-fast joining based on genomic coordinate and mutant base
setkey(obs_dt, position, new_base)
setkey(sim_dt, position, new_base)

# Perform an inner join to find all simulation instances where the observed mutation occurred
# merged_dt contains one row per simulation replicate for the specific mutated sites
merged_dt <- sim_dt[obs_dt, nomatch = 0]

# Count how many simulation replicates produced a homoplasy count >= the observed count
# (This fulfills the reviewer's strict "analyzed per replicate" requirement)
p_val_dt <- merged_dt[simulated_homoplasy_count >= count,
                      .(exceed_count = .N),
                      by = .(position, new_base)]

# Merge the exceed_count back into the original observed dataset
obs_dt <- merge(obs_dt, p_val_dt, by = c("position", "new_base"), all.x = TRUE)

# If a mutation was NEVER generated in 1000 simulations, exceed_count is NA. Convert to 0.
obs_dt[is.na(exceed_count), exceed_count := 0]

# Calculate the Empirical P-value using the standard pseudo-count method: (n + 1) / (N + 1)
# This avoids P=0, which is statistically inappropriate for Monte Carlo simulations.
obs_dt[, empirical_p_value := (exceed_count + 1) / (num_simulations + 1)]

# ------------------------------------------------------------------------------
# 5. Multiple Testing Correction (FDR - Benjamini-Hochberg)
# ------------------------------------------------------------------------------
cat("Applying Benjamini-Hochberg FDR correction...\n")
obs_dt[, fdr_adjusted_p := p.adjust(empirical_p_value, method = "BH")]

# Flag significantly selected variants (FDR < 0.05)
obs_dt[, is_significant := fdr_adjusted_p < 0.05]

# Sort the final table: Most significant first, then by highest observed homoplasy count
setorder(obs_dt, fdr_adjusted_p, empirical_p_value, -count)
obs_df <- as.data.frame(obs_dt) %>% select(-new_base) %>% 
    select(gene,position,codon_pos,ref,alt,type,count,exceed_count,empirical_p_value,fdr_adjusted_p,is_significant)
# ------------------------------------------------------------------------------
# 6. Save Final Results
# ------------------------------------------------------------------------------
output_file <- "final_convergent_mutations_statistics.csv"
fwrite(obs_df, output_file)

cat("Analysis complete! Statistically significant targets successfully identified.\n")
cat("Results saved to:", output_file, "\n")

# Preview the top significant hits
print(head(obs_dt[is_significant == TRUE]))