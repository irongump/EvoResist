# 02-make_initial_list.R
#
# Build the annotated initial candidate variant list for a single drug.
#
# Merges convergent SNP and indel calls from the 100k dataset with
# WHO catalogue variant names, assigns gene-body or promoter annotations,
# aggregates convergence event numbers by genomic group, and writes the
# combined list to disk.
#
# Usage (called from pipeline.sh):
#   Rscript ./02-make_initial_list.R <drug_name> <gene_list> \
#       <start_list> <end_list> <position_strand_list> <results_dir>
#
# Args (positional):
#   1  drug_name            e.g. RIF
#   2  gene_list            comma-separated gene names, e.g. inhA,katG
#   3  start_list           comma-separated gene start positions
#   4  end_list             comma-separated gene end positions
#   5  position_strand_list comma-separated strand flags (1=forward, 0=reverse)
#   6  results_dir          base results directory

args <- commandArgs(trailingOnly = TRUE)

drug_name            <- args[1]
gene_list            <- unlist(strsplit(args[2], ","))
start_list           <- as.integer(unlist(strsplit(args[3], ",")))
end_list             <- as.integer(unlist(strsplit(args[4], ",")))
position_strand_list <- as.integer(unlist(strsplit(args[5], ",")))
results_dir          <- args[6]

library(data.table)
library(dplyr)

# ------------------------------------------------------------------
# Load SNP data (100k convergent variants, filtered to this drug)
# ------------------------------------------------------------------
input_100k_snp <- fread(
  file.path(results_dir, drug_name, "denovo_snp_2.txt"),
  header = FALSE, fill = TRUE
)
setnames(input_100k_snp, c("gene", "codon", "anno", "position", "ref", "alt", "Event_number"))
input_100k_snp <- input_100k_snp[!(ref == "" | alt == "")]

# Load WHO catalogue variants for cross-referencing
input_who <- fread(
  file.path(results_dir, drug_name, "WHO_list", "WHO_list_allGroup.txt"),
  header = FALSE
)
setnames(input_who, c("variant_name_who", "position", "ref", "alt", "grading"))

gene_info <- data.table(
  target_gene = gene_list,
  start       = start_list,
  end         = end_list,
  strand      = position_strand_list
)

# ------------------------------------------------------------------
# Annotate variants with gene-body or promoter region labels.
# Multi-gene targets are handled via interval overlap; the nearest
# promoter gene is used for variants outside all coding regions.
# inhA has a special extended promoter (covers the fabG1 upstream region).
# rrs and rrl variants are always grouped by position.
# ------------------------------------------------------------------
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
    # Retain nearest promoter gene for each variant
    promoter_map <- promoter_map[, .SD[which.min(abs(rel_pos))], by = row_id][, .(row_id, promoter_name)]
    dt <- merge(dt, promoter_map, by = "row_id", all.x = TRUE, sort = FALSE)
  } else {
    dt[, promoter_name := NA_character_]
  }

  # Step 3: build final annotation name
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

# ------------------------------------------------------------------
# Aggregate convergence event numbers by genomic group.
# Coding variants are grouped by codon+annotation; promoter variants
# by position; rrs/rrl always by position.
# inhA has a hard-coded extended promoter boundary.
# ------------------------------------------------------------------
annotate_variant_group <- function(aggr_all_snp,
                                   gene_info,
                                   promoter_bp             = 1000,
                                   inhA_force_promoter     = TRUE,
                                   inhA_promoter_start     = 1672440,
                                   inhA_start              = 1674202,
                                   verbose                 = TRUE) {

  dt <- copy(as.data.table(aggr_all_snp))
  gi <- copy(as.data.table(gene_info))

  req_dt <- c("position", "ref", "alt", "Event_number", "gene", "codon", "anno")
  req_gi <- c("target_gene", "start", "end", "strand")

  if (length(setdiff(req_dt, names(dt))) > 0)
    stop("aggr_all_snp missing: ", paste(setdiff(req_dt, names(dt)), collapse = ", "))
  if (length(setdiff(req_gi, names(gi))) > 0)
    stop("gene_info missing: ", paste(setdiff(req_gi, names(gi)), collapse = ", "))

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

  # Deduplicate by position-ref-alt to avoid double-counting
  dt[, variant_key := paste(position, ref, alt, sep = ":")]
  dup_info <- dt[, .N, by = variant_key][N > 1]
  if (verbose && nrow(dup_info) > 0)
    message("Found ", nrow(dup_info), " duplicated position-ref-alt keys; keeping first occurrence only.")
  dt <- dt[!duplicated(variant_key)]

  # Parse single and range positions
  dt[, pos_start := as.numeric(sub("-.*$", "", position))]
  dt[, pos_end   := fifelse(grepl("-", position), as.numeric(sub("^.*-", "", position)), pos_start)]
  dt[, row_id__  := .I]

  dt[, `:=`(group_type = "other", matched_target_gene = NA_character_, group_id = NA_character_)]

  # rrs and rrl: always group by position
  dt[gene %chin% c("rrs", "rrl"),
     `:=`(group_type = "rrs_rrl", matched_target_gene = gene, group_id = paste(gene, pos_start, sep = "|"))]

  # Build promoter intervals
  prom_iv <- gi[, .(
    target_gene,
    start = fifelse(strand == 1, start - promoter_bp, end + 1),
    end   = fifelse(strand == 1, start - 1,           end + promoter_bp)
  )]
  prom_iv[, `:=`(start2 = pmin(start, end), end2 = pmax(start, end))]
  prom_iv[, `:=`(start = start2, end = end2)]
  prom_iv[, c("start2", "end2") := NULL]

  # Override inhA promoter to its known upstream boundary
  if (inhA_force_promoter && "inhA" %chin% prom_iv$target_gene) {
    prom_iv[target_gene == "inhA", `:=`(start = inhA_promoter_start, end = inhA_start - 1)]
  }

  coding_iv <- gi[, .(target_gene, start, end)]

  # Force-assign inhA promoter region first
  if (inhA_force_promoter && "inhA" %chin% prom_iv$target_gene) {
    inhA_iv <- prom_iv[target_gene == "inhA"]
    dt[group_type == "other" & pos_start >= inhA_iv$start[1] & pos_start <= inhA_iv$end[1],
       `:=`(group_type = "promoter", matched_target_gene = "inhA",
            group_id = paste("promoter", pos_start, sep = "|"))]
  }

  # Coding overlap
  x_coding <- dt[group_type == "other", .(row_id__, start = pos_start, end = pos_end, gene, codon, anno)]
  if (nrow(x_coding) > 0) {
    setkey(x_coding, start, end)
    setkey(coding_iv, start, end)
    coding_hit <- foverlaps(x_coding, coding_iv, by.x = c("start", "end"), by.y = c("start", "end"),
                            type = "any", nomatch = NULL)
    if (nrow(coding_hit) > 0) {
      coding_hit <- coding_hit[, .SD[1], by = row_id__]
      coding_hit[, `:=`(group_type = "coding", matched_target_gene = target_gene,
                         group_id = paste(gene, codon, anno, sep = "|"))]
      dt[coding_hit, on = "row_id__",
         `:=`(group_type = i.group_type, matched_target_gene = i.matched_target_gene, group_id = i.group_id)]
    }
  }

  # Promoter overlap for remaining unclassified variants
  x_prom <- dt[group_type == "other", .(row_id__, start = pos_start, end = pos_start)]
  if (nrow(x_prom) > 0) {
    setkey(x_prom, start, end)
    setkey(prom_iv, start, end)
    prom_hit <- foverlaps(x_prom, prom_iv, by.x = c("start", "end"), by.y = c("start", "end"),
                          type = "within", nomatch = NULL)
    if (nrow(prom_hit) > 0) {
      prom_hit2 <- prom_hit[, .(matched_target_gene = paste(unique(target_gene), collapse = ";")), by = row_id__]
      prom_hit2 <- dt[prom_hit2, on = "row_id__", .(row_id__, pos_start, matched_target_gene = i.matched_target_gene)]
      prom_hit2[, `:=`(group_type = "promoter", group_id = paste("promoter", pos_start, sep = "|"))]
      dt[prom_hit2, on = "row_id__",
         `:=`(group_type = i.group_type, matched_target_gene = i.matched_target_gene, group_id = i.group_id)]
    }
  }

  # Sum convergence event numbers within each group
  group_sum <- dt[!is.na(group_id), .(group_event_number = sum(Event_number, na.rm = TRUE)), by = group_id]
  dt[group_sum, on = "group_id", group_event_number := i.group_event_number]

  new_cols <- c("group_type", "matched_target_gene", "group_id", "group_event_number")
  out_cols <- c(names(aggr_all_snp), new_cols)
  dt <- dt[, ..out_cols]

  return(dt[])
}

# ------------------------------------------------------------------
# Separate continuous (MNP/range) and single-position WHO entries
# ------------------------------------------------------------------
who_continuous_snp <- input_who[grep("-", position), ]
who_single         <- input_who[!(grep("-", position)), ]
who_single[, position := as.integer(position)]

input_100k_continuous_snp <- input_100k_snp[grep("-", position), ]
input_100k_single_snp     <- input_100k_snp[!(grep("-", position)), ]
input_100k_single_snp[, position := as.integer(position)]

# Annotate single-position SNPs and merge with WHO names
input_100k_single_snp2 <- annotate_promoter_multi(input_100k_single_snp, gene_info)
combined_single_snp     <- merge(input_100k_single_snp2, who_single, by = c("position", "ref", "alt"), all.x = TRUE)
combined_single_snp[, variant := ifelse(is.na(variant_name_who), anno_name,
                                        paste0(variant_name_who, "(", anno_name, ")"))]
combined_single_snp_final <- combined_single_snp[!(ref == alt), ]

# Annotate continuous SNPs (MNPs / position ranges)
combined_continuous_snp  <- merge(input_100k_continuous_snp, who_continuous_snp, by = c("position", "ref", "alt"), all.x = TRUE)
combined_continuous_snp[, variant := ifelse(is.na(variant_name_who),
                                            paste(gene, codon, anno, sep = "_"),
                                            paste0(variant_name_who, "(", paste(gene, codon, anno, sep = "_"), ")"))]
combined_continuous_snp_final <- combined_continuous_snp[!(ref == alt), ]

output_snp_100k <- rbind(
  combined_single_snp_final[,    c("variant", "grading", "position", "ref", "alt", "Event_number")],
  combined_continuous_snp_final[, c("variant", "grading", "position", "ref", "alt", "Event_number")]
)

# ------------------------------------------------------------------
# Load and annotate indels
# ------------------------------------------------------------------
input_100k_indel <- fread(
  file.path(results_dir, drug_name, "denovo_indel_2.txt"),
  header = FALSE
)
setnames(input_100k_indel, c("Isolate_number", "position", "ref", "alt", "codon", "length_change", "gene"))
input_100k_indel <- input_100k_indel[!(ref == "" | alt == "")]

combined_indel <- merge(input_100k_indel, who_single, by = c("position", "ref", "alt"), all.x = TRUE)
combined_indel[, anno_name := paste(gene, codon, length_change, sep = "_")]
combined_indel[, variant   := ifelse(is.na(variant_name_who), anno_name,
                                     paste0(variant_name_who, "(", anno_name, ")"))]
combined_indel$Event_number <- NA
combined_indel <- combined_indel[!grepl(",", alt), ]
output_indel_100k <- combined_indel[, c("variant", "grading", "position", "ref", "alt", "Event_number")]

# ------------------------------------------------------------------
# Combine SNPs and indels, aggregate event numbers by variant group
# ------------------------------------------------------------------
output_list <- unique(rbind(output_snp_100k, output_indel_100k))
output_list[, grading      := ifelse(is.na(grading),      0, grading)]
output_list[, Event_number := ifelse(is.na(Event_number), 0, Event_number)]

aggr_all_snp <- rbind(
  combined_single_snp_final[,    c("variant", "grading", "position", "ref", "alt", "Event_number", "gene", "codon", "anno")],
  combined_continuous_snp_final[, c("variant", "grading", "position", "ref", "alt", "Event_number", "gene", "codon", "anno")]
)

aggregated_snp <- annotate_variant_group(
  aggr_all_snp = aggr_all_snp,
  gene_info    = gene_info,
  promoter_bp  = 1000
)

aggregated_all <- rbind(
  aggregated_snp,
  combined_indel[, c("variant", "grading", "position", "ref", "alt", "Event_number", "gene", "codon")],
  fill = TRUE
)

# ------------------------------------------------------------------
# Write annotated initial list
# ------------------------------------------------------------------
outdir <- file.path(results_dir, drug_name, "denovo_EvoResist_initial_list.txt")
fwrite(aggregated_all, outdir, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = NA)

message("Done: initial variant list written to ", outdir)
