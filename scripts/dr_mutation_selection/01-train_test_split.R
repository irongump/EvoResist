# 01-train_test_split.R
#
# Split a per-drug sample list into 70 % training and 30 % test sets,
# stratified by phenotype (R/S) so that resistance prevalence is
# preserved in both partitions.
#
# Usage (called from pipeline.sh):
#   Rscript ./01-train_test_split.R <sample_list_file> <output_dir>
#
# Arguments:
#   sample_list_file  Path to the tab-separated sample list, e.g.
#                     ./data/RIF_sample_list.txt
#                     Required columns: sample_id, pheno
#   output_dir        Directory where train_70.txt and test_30.txt
#                     will be written, e.g. <RESULTS_DIR>/RIF/id/
#
# Output files (tab-separated, with header):
#   <output_dir>/train_70.txt
#   <output_dir>/test_30.txt

library(data.table)

# ── Arguments ────────────────────────────────────────────────────────────────
args        <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 01-train_test_split.R <sample_list_file> <output_dir>")
}

sample_list_file <- args[1]
output_dir       <- args[2]

# ── Parameters ───────────────────────────────────────────────────────────────
train_ratio <- 0.7
set.seed(20260330)

# ── Load data ─────────────────────────────────────────────────────────────────
id_list <- fread(sample_list_file, header = TRUE, sep = "\t")

# ── Stratified split by phenotype ────────────────────────────────────────────
id_list[, row_id := .I]

train_rows <- id_list[, {
  n_train <- ceiling(.N * train_ratio)
  .(row_id = sample(row_id, n_train))
}, by = pheno]$row_id

train_dt <- id_list[ row_id %in% train_rows]
test_dt  <- id_list[!row_id %in% train_rows]

train_dt[, row_id := NULL]
test_dt[,  row_id := NULL]

# ── Write output ──────────────────────────────────────────────────────────────
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(train_dt, file.path(output_dir, "train_70.txt"),
       sep = "\t", quote = FALSE, na = NA, col.names = TRUE, row.names = FALSE)

fwrite(test_dt, file.path(output_dir, "test_30.txt"),
       sep = "\t", quote = FALSE, na = NA, col.names = TRUE, row.names = FALSE)

message("Done: ", nrow(train_dt), " train / ", nrow(test_dt), " test samples written to ", output_dir)
