# 02-make_variant_list.R
#
# Build the annotated initial candidate variant list for a single drug,
# then apply convergence-threshold and promoter-length filters to produce
# the filtered list1 for the given threshold × promoter combination.
#
# The initial list is built on the first call and cached on disk.
# Subsequent calls for the same drug (different threshold/promoter)
# skip the build step and read the cached file directly.
#
# Usage (called from pipeline.sh):
#   Rscript ./02-make_variant_list.R <drug_name> <gene_list> \
#       <start_list> <end_list> <LoF_list> <threshold> \
#       <promoter_region> <positive_strand_list> <results_dir>
#
# Args (positional):
#   1  drug_name             e.g. RIF
#   2  gene_list             comma-separated gene names, e.g. inhA,katG
#   3  start_list            comma-separated gene start positions
#   4  end_list              comma-separated gene end positions
#   5  LoF_list              comma-separated LoF flags (1=LoF gene, 0=not)
#   6  threshold             convergence event number threshold (integer)
#   7  promoter_region       promoter window in bp (integer)
#   8  positive_strand_list  comma-separated strand flags (1=forward, 0=reverse)
#   9  results_dir           base results directory

library(data.table)
library(dplyr)

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

gene_info <- data.table(
  target_gene = gene_list,
  start       = start_list,
  end         = end_list,
  strand      = positive_strand_list,
  LoF         = LoF_list
)

# ==================================================================
# PART 1 — Build annotated initial candidate list
# Skipped if the file already exists on disk (cache).
# ==================================================================

initial_list_path <- file.path(
  results_dir, drug_name, "denovo_EvoResist_initial_list.txt"
)

if (!file.exists(initial_list_path)) {

  message("Building initial list for ", drug_name, " ...")

  # ── Load SNP data ──────────────────────────────────────────────────────────
  input_100k_snp <- fread(
    file.path(results_dir, drug_name, "Sens_analysis_75", "denovo_snp_2.txt"),
    header = FALSE, fill = TRUE
  )
  setnames(input_100k_snp, c("gene", "codon", "anno", "position", "ref", "alt", "Event_number"))
  input_100k_snp <- input_100k_snp[!(ref == "" | alt == "")]

  # Load WHO catalogue for cross-referencing
  input_who <- fread(
    file.path(results_dir, drug_name, "WHO_list", "WHO_list_allGroup.txt"),
    header = FALSE
  )
  setnames(input_who, c("variant_name_who", "position", "ref", "alt", "grading"))

  # ── Promoter annotation ────────────────────────────────────────────────────
  annotate_promoter_multi <- function(dt, gene_info) {
    dt        <- copy(as.data.table(dt))
    gene_info <- copy(as.data.table(gene_info))

    dt[, row_id := .I]
    dt[, position := as.numeric(position)]
    gene_info[, `:=`(
      target_gene = as.character(target_gene),
      start       = as.numeric(start),
      end         = as.numeric(end),
      strand      = as.numeric(strand)
    )]

    # Step 1: identify variants inside a target gene body
    body_list <- lapply(seq_len(nrow(gene_info)), function(i) {
      g   <- gene_info[i]
      tmp <- dt[position >= g$start & position <= g$end, .(row_id, body_gene = g$target_gene)]
      if (nrow(tmp) == 0) return(NULL)
      tmp
    })
    body_map <- rbindlist(body_list, use.names = TRUE, fill = TRUE)

    if (nrow(body_map) > 0) {
      body_map <- body_map[, .(body_gene = paste(unique(body_gene), collapse = ";")), by = row_id]
      dt <- merge(dt, body_map, by = "row_id", all.x = TRUE, sort = FALSE)
    } else {
      dt[, body_gene := NA_character_]
    }

    # Step 2: assign promoter labels only to variants outside all gene bodies
    dt_nonbody <- dt[is.na(body_gene)]

    candidate_list <- lapply(seq_len(nrow(gene_info)), function(i) {
      g   <- gene_info[i]
      tmp <- copy(dt_nonbody)

      if (g$strand == 1) {
        tmp <- tmp[position < g$start]
        if (nrow(tmp) == 0) return(NULL)
        tmp[, rel_pos := position - g$start]
      } else if (g$strand == 0) {
        tmp <- tmp[position > g$end]
        if (nrow(tmp) == 0) return(NULL)
        tmp[, rel_pos := g$end - position]
      } else {
        stop("strand must contain only 1 or 0")
      }

      tmp[, promoter_name := paste0(g$target_gene, "_", rel_pos, ref, ">", alt)]
      tmp[, .(row_id, promoter_name, rel_pos)]
    })

    promoter_map <- rbindlist(candidate_list, use.names = TRUE, fill = TRUE)

    if (nrow(promoter_map) > 0) {
      promoter_map <- promoter_map[, .SD[which.min(abs(rel_pos))], by = row_id][, .(row_id, promoter_name)]
      dt <- merge(dt, promoter_map, by = "row_id", all.x = TRUE, sort = FALSE)
    } else {
      dt[, promoter_name := NA_character_]
    }

    dt[, anno_name := fifelse(
      !is.na(body_gene),
      paste(gene, codon, anno, sep = "_"),
      fifelse(
        !is.na(promoter_name),
        promoter_name,
        paste(gene, codon, anno, sep = "_")
      )
    )]

    dt[, c("row_id", "body_gene", "promoter_name") := NULL]
    return(dt[])
  }

  # ── Variant group aggregation ──────────────────────────────────────────────
  annotate_variant_group <- function(aggr_all_snp,
                                     gene_info,
                                     promoter_bp         = 1000,
                                     inhA_force_promoter = TRUE,
                                     inhA_promoter_start = 1672440,
                                     inhA_start          = 1674202,
                                     verbose             = TRUE) {

    dt <- copy(as.data.table(aggr_all_snp))
    gi <- copy(as.data.table(gene_info))

    dt[, `:=`(
      position     = as.character(position),
      ref          = as.character(ref),
      alt          = as.character(alt),
      gene         = as.character(gene),
      codon        = as.character(codon),
      anno         = as.character(anno),
      Event_number = as.numeric(Event_number)
    )]
    gi[, `:=`(
      target_gene = as.character(target_gene),
      start       = as.numeric(start),
      end         = as.numeric(end),
      strand      = as.numeric(strand)
    )]

    dt[, variant_key := paste(position, ref, alt, sep = ":")]
    dup_info <- dt[, .N, by = variant_key][N > 1]
    if (verbose && nrow(dup_info) > 0)
      message("Found ", nrow(dup_info), " duplicated position-ref-alt keys; keeping first occurrence only.")
    dt <- dt[!duplicated(variant_key)]

    dt[, pos_start := as.numeric(sub("-.*$", "", position))]
    dt[, pos_end   := fifelse(grepl("-", position), as.numeric(sub("^.*-", "", position)), pos_start)]
    dt[, row_id__  := .I]
    dt[, `:=`(group_type = "other", matched_target_gene = NA_character_, group_id = NA_character_)]

    dt[gene %chin% c("rrs", "rrl"),
       `:=`(group_type = "rrs_rrl", matched_target_gene = gene,
            group_id = paste(gene, pos_start, sep = "|"))]

    prom_iv <- gi[, .(
      target_gene,
      start = fifelse(strand == 1, start - promoter_bp, end + 1),
      end   = fifelse(strand == 1, start - 1,           end + promoter_bp)
    )]
    prom_iv[, `:=`(start2 = pmin(start, end), end2 = pmax(start, end))]
    prom_iv[, `:=`(start = start2, end = end2)][, c("start2", "end2") := NULL]

    if (inhA_force_promoter && "inhA" %chin% prom_iv$target_gene) {
      prom_iv[target_gene == "inhA", `:=`(start = inhA_promoter_start, end = inhA_start - 1)]
      inhA_iv <- prom_iv[target_gene == "inhA"]
      dt[group_type == "other" & pos_start >= inhA_iv$start[1] & pos_start <= inhA_iv$end[1],
         `:=`(group_type = "promoter", matched_target_gene = "inhA",
              group_id = paste("promoter", pos_start, sep = "|"))]
    }

    coding_iv <- gi[, .(target_gene, start, end)]
    x_coding  <- dt[group_type == "other", .(row_id__, start = pos_start, end = pos_end, gene, codon, anno)]
    if (nrow(x_coding) > 0) {
      setkey(x_coding, start, end); setkey(coding_iv, start, end)
      coding_hit <- foverlaps(x_coding, coding_iv, by.x = c("start","end"), by.y = c("start","end"),
                              type = "any", nomatch = NULL)
      if (nrow(coding_hit) > 0) {
        coding_hit <- coding_hit[, .SD[1], by = row_id__]
        coding_hit[, `:=`(group_type = "coding", matched_target_gene = target_gene,
                           group_id = paste(gene, codon, anno, sep = "|"))]
        dt[coding_hit, on = "row_id__",
           `:=`(group_type = i.group_type, matched_target_gene = i.matched_target_gene, group_id = i.group_id)]
      }
    }

    x_prom <- dt[group_type == "other", .(row_id__, start = pos_start, end = pos_start)]
    if (nrow(x_prom) > 0) {
      setkey(x_prom, start, end); setkey(prom_iv, start, end)
      prom_hit <- foverlaps(x_prom, prom_iv, by.x = c("start","end"), by.y = c("start","end"),
                            type = "within", nomatch = NULL)
      if (nrow(prom_hit) > 0) {
        prom_hit2 <- prom_hit[, .(matched_target_gene = paste(unique(target_gene), collapse = ";")), by = row_id__]
        prom_hit2 <- dt[prom_hit2, on = "row_id__",
                        .(row_id__, pos_start, matched_target_gene = i.matched_target_gene)]
        prom_hit2[, `:=`(group_type = "promoter", group_id = paste("promoter", pos_start, sep = "|"))]
        dt[prom_hit2, on = "row_id__",
           `:=`(group_type = i.group_type, matched_target_gene = i.matched_target_gene, group_id = i.group_id)]
      }
    }

    group_sum <- dt[!is.na(group_id), .(group_event_number = sum(Event_number, na.rm = TRUE)), by = group_id]
    dt[group_sum, on = "group_id", group_event_number := i.group_event_number]

    new_cols <- c("group_type", "matched_target_gene", "group_id", "group_event_number")
    out_cols <- c(names(aggr_all_snp), new_cols)
    return(dt[, ..out_cols])
  }

  # ── Merge SNPs with WHO names ──────────────────────────────────────────────
  who_continuous_snp        <- input_who[grep("-", position), ]
  who_single                <- input_who[!(grep("-", position)), ]
  who_single[, position     := as.integer(position)]

  input_100k_continuous_snp <- input_100k_snp[grep("-", position), ]
  input_100k_single_snp     <- input_100k_snp[!(grep("-", position)), ]
  input_100k_single_snp[, position := as.integer(position)]

  input_100k_single_snp2    <- annotate_promoter_multi(input_100k_single_snp, gene_info)
  combined_single_snp        <- merge(input_100k_single_snp2, who_single,
                                      by = c("position","ref","alt"), all.x = TRUE)
  combined_single_snp[, variant := ifelse(is.na(variant_name_who), anno_name,
                                          paste0(variant_name_who, "(", anno_name, ")"))]
  combined_single_snp_final  <- combined_single_snp[!(ref == alt), ]

  combined_continuous_snp    <- merge(input_100k_continuous_snp, who_continuous_snp,
                                      by = c("position","ref","alt"), all.x = TRUE)
  combined_continuous_snp[, variant := ifelse(is.na(variant_name_who),
                                              paste(gene, codon, anno, sep = "_"),
                                              paste0(variant_name_who, "(", paste(gene, codon, anno, sep = "_"), ")"))]
  combined_continuous_snp_final <- combined_continuous_snp[!(ref == alt), ]

  # ── Load and annotate indels ───────────────────────────────────────────────
  input_100k_indel <- fread(
    file.path(results_dir, drug_name, "Sens_analysis_75", "denovo_indel_2.txt"),
    header = FALSE
  )
  setnames(input_100k_indel, c("Isolate_number","position","ref","alt","codon","length_change","gene"))
  input_100k_indel <- input_100k_indel[!(ref == "" | alt == "")]

  combined_indel <- merge(input_100k_indel, who_single,
                          by = c("position","ref","alt"), all.x = TRUE)
  combined_indel[, anno_name := paste(gene, codon, length_change, sep = "_")]
  combined_indel[, variant   := ifelse(is.na(variant_name_who), anno_name,
                                       paste0(variant_name_who, "(", anno_name, ")"))]
  combined_indel$Event_number <- NA
  combined_indel <- combined_indel[!grepl(",", alt), ]

  # ── Aggregate event numbers by variant group ───────────────────────────────
  aggr_all_snp <- rbind(
    combined_single_snp_final[,    c("variant","grading","position","ref","alt","Event_number","gene","codon","anno")],
    combined_continuous_snp_final[, c("variant","grading","position","ref","alt","Event_number","gene","codon","anno")]
  )

  aggregated_snp <- annotate_variant_group(aggr_all_snp, gene_info, promoter_bp = 1000)

  aggregated_all <- rbind(
    aggregated_snp,
    combined_indel[, c("variant","grading","position","ref","alt","Event_number","gene","codon")],
    fill = TRUE
  )
  aggregated_all[, grading      := ifelse(is.na(grading),      0, grading)]
  aggregated_all[, Event_number := ifelse(is.na(Event_number), 0, Event_number)]

  dir.create(dirname(initial_list_path), recursive = TRUE, showWarnings = FALSE)
  fwrite(aggregated_all, initial_list_path,
         col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = NA)
  message("Initial list written to ", initial_list_path)

} else {
  message("Initial list already exists for ", drug_name, ", skipping rebuild.")
}

# ==================================================================
# PART 2 — Apply threshold × promoter filter to produce list1
# ==================================================================

message("Applying threshold=", threshold, ", promoter=", promoter_region, " for ", drug_name, " ...")

initial_list  <- fread(initial_list_path)
initial_snp   <- initial_list[!is.na(Event_number) & Event_number != 0, ]
initial_indel <- initial_list[ is.na(Event_number) | Event_number == 0, ]

print(paste0("Drug: ", drug_name, ". SNPs: ", nrow(initial_snp), ", Indels: ", nrow(initial_indel)))

# ── Indel selection ────────────────────────────────────────────────────────
# LoF genes: include all indels (frameshift and in-frame).
# Non-LoF genes: include in-frame indels only.
# rrs / rrl: no indels included.
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

# ── SNP selection ──────────────────────────────────────────────────────────
# Promoter variants: within the promoter window with event number >= threshold.
# Coding variants: non-synonymous with event number >= threshold;
#                  premature-stop added unconditionally for LoF genes.
# inhA: extended promoter covering the fabG1 upstream region.
# gyrA: promoter restricted to the intergenic region between gyrB and gyrA.
# STM: two specific rrs positions always evaluated.
initial_snp_2 <- copy(initial_snp)
initial_snp_2[, pos_start := as.numeric(sub("-.*$", "", position))]
initial_snp_2[, pos_end   := fcase(grepl("-", position),
                                    as.numeric(sub("^.*-", "", position)),
                                    default = pos_start)]

snps_incldued_list <- list()
for (i in seq_along(gene_list)) {
  gene_name            <- gene_list[i]
  gene_start           <- start_list[i]
  gene_end             <- end_list[i]
  gene_LoF             <- (LoF_list[i] == 1)
  gene_positive_strand <- (positive_strand_list[i] == 1)

  if (gene_name == "inhA") {
    inhA_promoter_start    <- gene_start - promoter_region
    SNPs_promoter_inhA     <- initial_snp_2[(pos_start >= inhA_promoter_start) & (pos_start < gene_start), ]
    fabG1_start            <- 1673440
    fabG1_promoter_start   <- fabG1_start - promoter_region
    SNPs_promoter_fabG1    <- initial_snp_2[(pos_start >= fabG1_promoter_start) & (pos_start < fabG1_start), ]
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

  SNPs_body_thisgene       <- initial_snp_2[(pos_start >= gene_start) & (pos_start <= gene_end), ]
  SNPs_body_thisgene_nosyn <- SNPs_body_thisgene[!(grepl("^Synonymous-", anno)), ]
  prestop_thisgene         <- SNPs_body_thisgene_nosyn[grepl("-stop$", anno), ]
  SNPs_event_pass_thisgene <- SNPs_body_thisgene_nosyn[group_event_number >= threshold, ]

  if (gene_LoF) {
    SNPs_body_thisgene_keep <- unique(rbind(SNPs_event_pass_thisgene, prestop_thisgene))
  } else {
    SNPs_body_thisgene_keep <- SNPs_event_pass_thisgene[!grepl("-stop$", anno), ]
  }

  snps_incldued_list[[gene_name]] <- unique(rbind(SNPs_body_thisgene_keep, SNPs_promoter_thisgene_keep))
}
snps_incldued <- do.call(rbind, snps_incldued_list)

# STM: two specific rrs positions always evaluated
if (drug_name == "STM") {
  STM_rrs       <- initial_snp_2[(pos_start %in% c(1472359, 1472362)) & (group_event_number >= threshold), ]
  snps_incldued <- rbind(snps_incldued, STM_rrs)
}

# ── Write list1 ────────────────────────────────────────────────────────────
output_list_cols <- c("variant","grading","position","ref","alt","Event_number",
                      "gene","codon","group_type","matched_target_gene","group_id","group_event_number")

out_list <- rbind(
  snps_incldued[,   ..output_list_cols],
  indels_incldued[, ..output_list_cols]
)

list1_path <- file.path(
  results_dir, drug_name, 
  paste0("Threshold_", threshold, "_Promoter_", promoter_region, "_list1.tsv")
)

fwrite(out_list, list1_path,
       col.names = TRUE, row.names = FALSE, na = NA, sep = "\t", quote = FALSE)

message("list1 written to ", list1_path)
