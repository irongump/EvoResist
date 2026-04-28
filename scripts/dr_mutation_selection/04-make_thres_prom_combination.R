# 04-make_thres_prom_combination.R
#
# Generate variant list1 for a given threshold × promoter-length combination.
#
# Selects variants from the annotated initial candidate list according to:
#   SNPs : convergence event number >= threshold for promoter/coding variants;
#          premature-stop mutations included unconditionally for LoF genes.
#   Indels: frameshift indels included for LoF genes; in-frame indels for
#           non-LoF genes; rrs/rrl excluded entirely.
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

# Load the annotated initial candidate list
initial_list  <- fread(file.path(results_dir, drug_name, "denovo_EvoResist_initial_list.txt"))
initial_snp   <- initial_list[!is.na(Event_number), ]
initial_indel <- initial_list[is.na(Event_number),  ]

print(paste0("Drug: ", drug_name, ". SNPs: ", nrow(initial_snp), ", Indels: ", nrow(initial_indel),
             ". nrow matching: ", (nrow(initial_snp) + nrow(initial_indel) == nrow(initial_list))))

# ------------------------------------------------------------------
# Indel selection
# LoF genes: include frameshift indels directly; evaluate in-frame below.
# Non-LoF genes: evaluate in-frame indels only.
# rrs / rrl: no indels included.
# ------------------------------------------------------------------
initial_indel[, fs       := abs(nchar(ref) - nchar(alt)) %% 3 != 0]
initial_indel[, position := as.integer(position)]

indels_incldued_list <- list()
for (i in seq_along(gene_list)) {
  gene_name  <- gene_list[i]
  gene_start <- start_list[i]
  gene_end   <- end_list[i]
  gene_LoF   <- (LoF_list[i] == 1)

  if (gene_name %in% c("rrs", "rrl")) {
    indels_incldued_list[[gene_name]] <- initial_indel[0]
    next
  }

  indel_thisgene <- initial_indel[(position >= gene_start) & (position <= gene_end), ]

  if (gene_LoF) {
    indels_incldued_list[[gene_name]] <- indel_thisgene
  } else {
    indels_incldued_list[[gene_name]] <- indel_thisgene[!(fs)]
  }
}
indels_incldued <- rbindlist(indels_incldued_list, use.names = TRUE, fill = TRUE)

# ------------------------------------------------------------------
# SNP selection
# Promoter: variants within [gene_start - promoter_region, gene_start)
#           (or symmetric on the reverse strand) with event number >= threshold.
# Coding:   non-synonymous variants with event number >= threshold;
#           premature-stop added unconditionally for LoF genes.
# inhA has a special two-boundary promoter (inhA + fabG1 upstream).
# gyrA promoter is restricted to the intergenic region between gyrB and gyrA.
# STM has two specific rrs positions that are always evaluated.
# ------------------------------------------------------------------
initial_snp_2 <- initial_snp[, pos_start := as.numeric(sub("-.*$", "", position))]
initial_snp_2[, pos_end := fcase(grepl("-", position), as.numeric(sub("^.*-", "", position)), default = pos_start)]

snps_incldued_list <- list()
for (i in seq_along(gene_list)) {
  gene_name            <- gene_list[i]
  gene_start           <- start_list[i]
  gene_end             <- end_list[i]
  gene_LoF             <- (LoF_list[i] == 1)
  gene_positive_strand <- (positive_strand_list[i] == 1)

  if (gene_name == "inhA") {
    inhA_promoter_start  <- gene_start - promoter_region
    SNPs_promoter_inhA   <- initial_snp_2[(pos_start >= inhA_promoter_start) & (pos_start < gene_start), ]
    fabG1_start          <- 1673440
    fabG1_promoter_start <- fabG1_start - promoter_region
    SNPs_promoter_fabG1  <- initial_snp_2[(pos_start >= fabG1_promoter_start) & (pos_start < fabG1_start), ]
    SNPs_promoter_thisgene <- unique(rbind(SNPs_promoter_inhA, SNPs_promoter_fabG1, fill = TRUE))
  } else if (gene_name == "gyrA") {
    SNPs_promoter_thisgene <- initial_snp_2[pos_start >= 7268 & pos_start <= 7301, ]
  } else if (gene_positive_strand) {
    promoter_start_thisgene <- gene_start - promoter_region
    SNPs_promoter_thisgene  <- initial_snp_2[(pos_start >= promoter_start_thisgene) & (pos_start < gene_start), ]
  } else {
    promoter_start_thisgene <- gene_end + promoter_region
    SNPs_promoter_thisgene  <- initial_snp_2[(pos_start <= promoter_start_thisgene) & (pos_start > gene_end), ]
  }

  SNPs_promoter_thisgene_keep <- SNPs_promoter_thisgene[group_event_number >= threshold, ]

  SNPs_body_thisgene        <- initial_snp_2[(pos_start >= gene_start) & (pos_start <= gene_end), ]
  SNPs_body_thisgene_nosyn  <- SNPs_body_thisgene[!(grepl("^Synonymous-", anno)), ]
  prestop_thisgene          <- SNPs_body_thisgene_nosyn[grepl("-stop$", anno), ]
  SNPs_event_pass_thisgene  <- SNPs_body_thisgene_nosyn[group_event_number >= threshold, ]

  if (gene_LoF) {
    SNPs_body_thisgene_keep <- unique(rbind(SNPs_event_pass_thisgene, prestop_thisgene))
  } else {
    SNPs_body_thisgene_keep <- SNPs_event_pass_thisgene[!grepl("-stop$", anno), ]
  }

  snps_incldued_list[[gene_name]] <- unique(rbind(SNPs_body_thisgene_keep, SNPs_promoter_thisgene_keep))
}
snps_incldued <- do.call(rbind, snps_incldued_list)

# STM: two specific rrs positions evaluated separately
if (drug_name == "STM") {
  STM_rrs       <- initial_snp_2[(pos_start %in% c(1472359, 1472362)) & (group_event_number >= threshold), ]
  snps_incldued <- rbind(snps_incldued, STM_rrs)
}

# ------------------------------------------------------------------
# Write list1
# ------------------------------------------------------------------
output_list_cols <- c("variant", "grading", "position", "ref", "alt", "Event_number",
                      "gene", "codon", "group_type", "matched_target_gene", "group_id", "group_event_number")

out_list <- rbind(
  snps_incldued[,    ..output_list_cols],
  indels_incldued[,  ..output_list_cols]
)

fwrite(
  out_list,
  file.path(results_dir, drug_name,
            paste0("Threshold_", threshold, "_Promoter_", promoter_region, "_list1.tsv")),
  col.names = TRUE, row.names = FALSE, na = NA, sep = "\t", quote = FALSE
)
