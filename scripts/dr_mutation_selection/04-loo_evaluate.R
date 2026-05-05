# 05-loo_evaluate.R
#
# Remove variants from list1 whose leave-one-out performance on the
# training set shows specificity harm exceeding sensitivity gain.
#
# Reads the per-variant leave-one-out analysis produced by
# 03-leave_one_out.py on the training set (train_list1/), then
# applies the following selection rules to produce list2:
#
#   Frameshift indels in LoF genes  -> retained unconditionally.
#   Premature-stop SNPs in LoF genes -> retained unconditionally.
#   All other variants               -> retained only if removing them
#                                       does not improve the combined
#                                       sensitivity + specificity score
#                                       (balance_weighted_contribution >= 0).
#
# Expert overrides (drug-specific):
#   RIF: rpoB_p.Leu430Pro (761095:T:C) always included.
#   LFX/MFX: gyrA-34C>T (7268:C:T) and gyrB_-165C>T (5075:C:T) always excluded.
#
# Args (positional):
#   1  drug_name             e.g. RIF
#   2  gene_list             comma-separated gene names
#   3  start_list            comma-separated gene start positions
#   4  end_list              comma-separated gene end positions
#   5  LoF_list              comma-separated LoF flags (1=LoF gene, 0=not)
#   6  threshold             convergence event number threshold (integer)
#   7  promoter_region       promoter window in bp (integer)
#   8  positive_strand_list  comma-separated strand flags (1=forward, 0=reverse)
#   9  results_dir           base results directory

library(data.table)

args <- commandArgs(trailingOnly = TRUE)

drug_name            <- args[1]
gene_list            <- unlist(strsplit(args[2], ","))
start_list           <- as.integer(unlist(strsplit(args[3], ",")))
end_list             <- as.integer(unlist(strsplit(args[4], ",")))
LoF_list             <- as.integer(unlist(strsplit(args[5], ",")))
threshold            <- as.integer(args[6])
promoter_region      <- as.integer(args[7])
positive_strand_list <- as.integer(unlist(strsplit(args[8], ",")))
results_dir          <- args[9]

print(paste0("Drug: ", drug_name, ", threshold ", threshold, ", promoter ", promoter_region, "."))

# Load per-variant leave-one-out results from the training set evaluation
combo_dir    <- file.path(results_dir, drug_name,
                          paste0("Threshold_", threshold, "_Promoter_", promoter_region))
list1_result <- fread(file.path(combo_dir, "train_list1", "per_variant_analysis.tsv"))

# Balance-weighted contribution: positive means the variant helps overall
# (sensitivity + specificity both increase when it is present)
list1_result[, balance_weighted_contribution := -(Delta_Sensitivity + Delta_Specificity)]

# Separate SNPs (have Event_number) from indels (Event_number is NA)
presensana_list_snp   <- list1_result[!is.na(Event_number), ]
presensana_list_indel <- list1_result[is.na(Event_number),  ]

# ------------------------------------------------------------------
# Indel selection
# ------------------------------------------------------------------
presensana_list_indel[, fs := abs(nchar(Ref) - nchar(Alt)) %% 3 != 0]

indels_LoF_list       <- list()
indels_selection_list <- list()

for (i in seq_along(gene_list)) {
  gene_name  <- gene_list[i]
  gene_start <- start_list[i]
  gene_end   <- end_list[i]
  gene_LoF   <- (LoF_list[i] == 1)

  indel_thisgene <- presensana_list_indel[(Position >= gene_start) & (Position <= gene_end), ]

  if (gene_LoF) {
    indels_LoF_list[[gene_name]]       <- indel_thisgene[(fs)]
    indels_selection_list[[gene_name]] <- indel_thisgene[!(fs)]
  } else {
    indels_selection_list[[gene_name]] <- indel_thisgene[!(fs)]
  }
}

indels_LoF       <- do.call(rbind, indels_LoF_list)
indels_selection <- do.call(rbind, indels_selection_list)

indels_selection_pass     <- indels_selection[!(balance_weighted_contribution < 0), ]
indels_selection_NOTpass  <- indels_selection[(balance_weighted_contribution < 0),  ]

# ------------------------------------------------------------------
# SNP selection
# ------------------------------------------------------------------
presensana_list_snp[, prestop := grepl("-stop$", group_id)]
presensana_list_snp_2 <- presensana_list_snp[, pos_start := as.numeric(sub("-.*$", "", Position))]
presensana_list_snp_2[, pos_end := fcase(grepl("-", Position), as.numeric(sub("^.*-", "", Position)), default = pos_start)]

SNPs_LoF_list       <- list()
SNPs_selection_list <- list()

for (i in seq_along(gene_list)) {
  gene_name  <- gene_list[i]
  gene_LoF   <- (LoF_list[i] == 1)

  snp_thisgene <- presensana_list_snp_2[matched_target_gene == gene_name, ]

  if (gene_LoF) {
    SNPs_LoF_list[[gene_name]]       <- snp_thisgene[(prestop)]
    SNPs_selection_list[[gene_name]] <- snp_thisgene[!(prestop)]
  } else {
    SNPs_selection_list[[gene_name]] <- snp_thisgene[!(prestop)]
  }
}

SNPs_LoF       <- do.call(rbind, SNPs_LoF_list)
SNPs_selection <- do.call(rbind, SNPs_selection_list)

# STM: two specific rrs positions evaluated alongside the main selection
if (drug_name == "STM") {
  STM_rrs        <- presensana_list_snp_2[(pos_start %in% c(1472359, 1472362)) & (group_event_number >= threshold), ]
  SNPs_selection <- rbind(SNPs_selection, STM_rrs)
}

SNPs_selection_pass    <- SNPs_selection[!(balance_weighted_contribution < 0), ]
SNPs_selection_NOTpass <- SNPs_selection[(balance_weighted_contribution < 0),  ]

# ------------------------------------------------------------------
# Combine retained variants into list2
# ------------------------------------------------------------------
output_list_cols <- c("Variant", "grading", "Position", "Ref", "Alt", "Event_number",
                      "gene", "codon", "group_type", "matched_target_gene", "group_id", "group_event_number")

safe_select <- function(x, cols) {
  if (is.null(x) || nrow(x) == 0) return(data.table())
  x[, ..cols]
}

final_list <- rbindlist(list(
  safe_select(SNPs_selection_pass,    output_list_cols),
  safe_select(SNPs_LoF,               output_list_cols),
  safe_select(indels_selection_pass,  output_list_cols),
  safe_select(indels_LoF,             output_list_cols)
), use.names = TRUE, fill = TRUE)

# Expert overrides
if (drug_name == "RIF") {
  extra <- list1_result[Position == 761095 & Ref == "T" & Alt == "C", ..output_list_cols]
  final_list <- rbind(final_list, extra)
}
if (drug_name %in% c("LFX", "MFX")) {
  final_list <- final_list[!((Position == 7268 & Ref == "C" & Alt == "T") |
                              (Position == 5075 & Ref == "C" & Alt == "T")), ]
}

# ------------------------------------------------------------------
# Write list2 and removed-variant record
# ------------------------------------------------------------------
setnames(final_list, c("variant", "grading", "position", "ref", "alt", "Event_number",
                        "gene", "codon", "group_type", "matched_target_gene", "group_id", "group_event_number"))

fwrite(
  final_list,
  file.path(results_dir, drug_name,
            paste0("Threshold_", threshold, "_Promoter_", promoter_region, "_list2.tsv")),
  col.names = TRUE, sep = "\t", quote = FALSE, na = NA
)
