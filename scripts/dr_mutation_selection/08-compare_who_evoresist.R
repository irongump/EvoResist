# 08-compare_who_evoresist.R
#
# Compare diagnostic performance of three variant catalogues:
#   G1       WHO Grade 1 only
#   G12      WHO Grade 1 + Grade 2
#   EvoResist  EvoResist final list (best threshold × promoter)
#
# For each drug, reads isolate-level predictions produced by
# 03-leave_one_out.py from the Final_Evaluate directories, then
# computes sensitivity, specificity, PPV, and NPV for each catalogue
# and their pairwise deltas (EvoResist vs G1 and EvoResist vs G12)
# with 95% CIs (Wilson for point estimates; Newcombe paired + bootstrap
# for deltas of sens/spec; bootstrap for PPV/NPV deltas).
# McNemar p-values are BH-adjusted across all 13 drugs.
#
# Usage (called from pipeline.sh):
#   Rscript ./08-compare_who_evoresist.R <results_dir> <output_dir> [group_col]
#
# Arguments:
#   results_dir   base results directory (same as RESULTS_DIR in pipeline.sh)
#   output_dir    directory where output TSV files are written
#   group_col     optional: column in isolate_predictions.tsv used for
#                 stratified analysis (e.g. Lineage). Must be present in all
#                 three prediction files. When supplied the script runs both
#                 an overall and a per-stratum analysis.
#
# Output files (always written):
#   <output_dir>/comparison_formatted.tsv
# Additional files (only when group_col is supplied):
#   <output_dir>/comparison_stratified_formatted.tsv

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 08-compare_who_evoresist.R <results_dir> <output_dir> [group_col]")
}
results_dir <- args[1]
output_dir  <- args[2]
group_col   <- if (length(args) >= 3) args[3] else NULL   # optional stratification column
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

drug_list <- c("RIF","INH","EMB","PZA","LFX","MFX","BDQ","LZD","AMK","STM","ETO","KAN","CAP")

# ── CI helpers ────────────────────────────────────────────────────────────────

wilson_ci <- function(x, n, conf.level = 0.95) {
  if (is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  z      <- qnorm(1 - (1 - conf.level) / 2)
  phat   <- x / n
  denom  <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  half   <- (z * sqrt((phat * (1 - phat) + z^2 / (4 * n)) / n)) / denom
  c(max(0, center - half), min(1, center + half))
}

newcombe_paired_ci <- function(n01, n10, N, conf.level = 0.95) {
  ci01 <- wilson_ci(n01, N, conf.level)
  ci10 <- wilson_ci(n10, N, conf.level)
  c(ci01[1] - ci10[2], ci01[2] - ci10[1])
}

bootstrap_paired_ci <- function(dt_sub, pred_a, pred_b,
                                conf.level = 0.95, B = 5000, seed = 1,
                                correct_label = "R") {
  set.seed(seed)
  N <- nrow(dt_sub)
  if (N == 0) return(c(NA_real_, NA_real_))
  a      <- dt_sub[[pred_a]]
  b      <- dt_sub[[pred_b]]
  deltas <- vapply(seq_len(B), function(i) {
    s   <- sample(N, N, replace = TRUE)
    n01 <- sum(a[s] != correct_label & b[s] == correct_label, na.rm = TRUE)
    n10 <- sum(a[s] == correct_label & b[s] != correct_label, na.rm = TRUE)
    (n01 - n10) / N
  }, numeric(1))
  alpha <- 1 - conf.level
  quantile(deltas, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE)
}

bootstrap_metric_diff <- function(dt, pred_a, pred_b, metric,
                                  conf.level = 0.95, B = 5000, seed = 1,
                                  positive_label = "R", pheno_col = "phenotype") {
  set.seed(seed)
  N <- nrow(dt)
  if (N == 0) return(c(NA_real_, NA_real_))
  deltas <- vapply(seq_len(B), function(i) {
    dt_s  <- dt[sample(N, N, replace = TRUE)]
    a_est <- calc_metric_single(dt_s, pred_a, metric,
                                positive_label = positive_label, pheno_col = pheno_col)$estimate
    b_est <- calc_metric_single(dt_s, pred_b, metric,
                                positive_label = positive_label, pheno_col = pheno_col)$estimate
    b_est - a_est
  }, numeric(1))
  alpha <- 1 - conf.level
  quantile(deltas, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE)
}

# ── Metric calculation ────────────────────────────────────────────────────────

calc_metric_single <- function(dt, pred_col, metric,
                               conf.level = 0.95,
                               positive_label = "R",
                               pheno_col = "phenotype") {
  pred  <- dt[[pred_col]]
  truth <- dt[[pheno_col]]
  denom_idx <- switch(metric,
                      sens = truth == positive_label,
                      spec = truth != positive_label,
                      ppv  = pred  == positive_label,
                      npv  = pred  != positive_label,
                      stop("Unknown metric: ", metric)
  )
  x   <- switch(metric,
                sens = sum(pred[denom_idx]  == positive_label, na.rm = TRUE),
                spec = sum(pred[denom_idx]  != positive_label, na.rm = TRUE),
                ppv  = sum(truth[denom_idx] == positive_label, na.rm = TRUE),
                npv  = sum(truth[denom_idx] != positive_label, na.rm = TRUE)
  )
  n   <- sum(denom_idx, na.rm = TRUE)
  est <- if (n == 0) NA_real_ else x / n
  list(x = x, n = n, estimate = est, ci = wilson_ci(x, n, conf.level))
}

compare_rate_paired <- function(dt, pred_a, pred_b,
                                truth_label, correct_label,
                                conf.level = 0.95, B = 5000, seed = 1,
                                alternative = "two.sided",
                                pheno_col = "phenotype") {
  dt_sub <- dt[get(pheno_col) == truth_label]
  N <- nrow(dt_sub)
  if (N == 0) return(list(N=0, n01=NA_integer_, n10=NA_integer_,
                          delta=NA_real_, ci_new=c(NA,NA), ci_boot=c(NA,NA), p=NA_real_))
  n01   <- sum(dt_sub[[pred_a]] != correct_label & dt_sub[[pred_b]] == correct_label, na.rm = TRUE)
  n10   <- sum(dt_sub[[pred_a]] == correct_label & dt_sub[[pred_b]] != correct_label, na.rm = TRUE)
  delta <- (n01 - n10) / N
  pval  <- if ((n01 + n10) == 0) NA_real_ else
    binom.test(n01, n01 + n10, alternative = alternative)$p.value
  list(N       = N,
       n01     = n01,
       n10     = n10,
       delta   = delta,
       ci_new  = newcombe_paired_ci(n01, n10, N, conf.level),
       ci_boot = bootstrap_paired_ci(dt_sub, pred_a, pred_b,
                                     conf.level, B, seed, correct_label),
       p       = pval)
}

compare_metric <- function(dt, pred_a, pred_b, metric,
                           conf.level = 0.95, B = 5000, seed = 1,
                           alternative = "two.sided",
                           positive_label = "R", negative_label = "S",
                           pheno_col = "phenotype") {
  a <- calc_metric_single(dt, pred_a, metric, conf.level, positive_label, pheno_col)
  b <- calc_metric_single(dt, pred_b, metric, conf.level, positive_label, pheno_col)
  
  if (metric %in% c("sens", "spec")) {
    truth_label   <- if (metric == "sens") positive_label else negative_label
    correct_label <- if (metric == "sens") positive_label else negative_label
    paired <- compare_rate_paired(dt, pred_a, pred_b, truth_label, correct_label,
                                  conf.level, B, seed, alternative, pheno_col)
    return(c(list(a=a, b=b), paired))
  }
  list(a=a, b=b, N=NA_integer_, n01=NA_integer_, n10=NA_integer_,
       delta   = b$estimate - a$estimate,
       ci_new  = c(NA_real_, NA_real_),
       ci_boot = bootstrap_metric_diff(dt, pred_a, pred_b, metric,
                                       conf.level, B, seed, positive_label, pheno_col),
       p       = NA_real_)
}

# ── Merge prediction tables from three catalogues ─────────────────────────────

standardize_pred_table <- function(dt, method_name,
                                   id_col, pheno_col, pred_col, drug_col,
                                   group_col = NULL) {
  dt   <- as.data.table(copy(dt))
  keep <- intersect(c(id_col, pheno_col, pred_col, drug_col, group_col), names(dt))
  dt   <- dt[, ..keep]
  setnames(dt, c(id_col, pheno_col, pred_col), c("id", "phenotype", paste0(method_name, "_pred")))
  if (!is.null(drug_col) && drug_col %in% names(dt)) setnames(dt, drug_col, "Drug")
  dt
}

merge_prediction_tables <- function(whoG1, G12, EvoResist,
                                    id_col, pheno_col, pred_col, drug_col,
                                    group_col = NULL) {
  g1   <- standardize_pred_table(whoG1,     "G1",   id_col, pheno_col, pred_col, drug_col, group_col)
  g12  <- standardize_pred_table(G12,       "G12",  id_col, pheno_col, pred_col, drug_col, group_col)
  evor <- standardize_pred_table(EvoResist, "EvoR", id_col, pheno_col, pred_col, drug_col, group_col)
  
  # group_col is retained from all three tables and added to merge_cols so it
  # passes through the join without duplication; all three files share the same
  # isolates so the group values are identical.
  merge_cols <- "id"
  if ("Drug" %in% names(g1) & "Drug" %in% names(g12) & "Drug" %in% names(evor))
    merge_cols <- c(merge_cols, "Drug")
  if (!is.null(group_col) && group_col %in% names(evor))
    merge_cols <- c(merge_cols, group_col)
  
  dt <- merge(g1, g12,   by = c(merge_cols, "phenotype"), all = FALSE)
  dt <- merge(dt,  evor, by = c(merge_cols, "phenotype"), all = FALSE)
  dt
}

# ── Build one result row per drug ─────────────────────────────────────────────

append_metric_block <- function(row, prefix, cmp_g1, cmp_g12) {
  for (method in c("G1","G12","EvoR")) {
    src <- if (method == "EvoR") cmp_g1$b else if (method == "G1") cmp_g1$a else cmp_g12$a
    row[, paste0(prefix,"_",method)          := src$estimate]
    row[, paste0(prefix,"_",method,"_CI_low"):= src$ci[1]]
    row[, paste0(prefix,"_",method,"_CI_high"):= src$ci[2]]
  }
  for (base in c("G1","G12")) {
    cmp <- if (base == "G1") cmp_g1 else cmp_g12
    row[, paste0(prefix,"_Delta_",base)              := cmp$delta]
    row[, paste0(prefix,"_Delta_",base,"_CI_Boot_low") := cmp$ci_boot[1]]
    row[, paste0(prefix,"_Delta_",base,"_CI_Boot_high"):= cmp$ci_boot[2]]
    if (prefix %in% c("Sens","Spec")) {
      row[, paste0(prefix,"_Delta_",base,"_CI_Newcombe_low") := cmp$ci_new[1]]
      row[, paste0(prefix,"_Delta_",base,"_CI_Newcombe_high"):= cmp$ci_new[2]]
      row[, paste0(prefix,"_P_",base)                        := cmp$p]
    }
  }
  row
}

build_result_row <- function(dt_sub, drug_name,
                             conf.level=0.95, B=5000, seed=20260220,
                             alternative="two.sided",
                             positive_label="R", negative_label="S") {
  row <- data.table(Drug = drug_name,
                    N_total = nrow(dt_sub),
                    N_R = sum(dt_sub$phenotype == positive_label, na.rm = TRUE))
  for (metric in c("sens","spec","ppv","npv")) {
    prefix <- toupper(substr(metric, 1, 1))
    prefix <- switch(metric, sens="Sens", spec="Spec", ppv="PPV", npv="NPV")
    cmp_g1  <- compare_metric(dt_sub, "G1_pred",  "EvoR_pred", metric,
                              conf.level, B, seed, alternative,
                              positive_label, negative_label)
    cmp_g12 <- compare_metric(dt_sub, "G12_pred", "EvoR_pred", metric,
                              conf.level, B, seed, alternative,
                              positive_label, negative_label)
    row <- append_metric_block(row, prefix, cmp_g1, cmp_g12)
  }
  row
}

# ── Formatting helpers ────────────────────────────────────────────────────────

fmt_pct_ci <- function(x, lo, hi, d=1) {
  ifelse(is.na(x), NA_character_,
         sprintf(paste0("%.",d,"f (%.",d,"f, %.",d,"f)"), 100*x, 100*lo, 100*hi))
}

fmt_p <- function(x, d=1) ifelse(is.na(x), NA_character_, sprintf(paste0("%.",d,"e"), x))

safe_col <- function(dt, col) if (col %in% names(dt)) dt[[col]] else rep(NA_real_, nrow(dt))

format_result_table <- function(dt, group_col=NULL, d=1, p_d=1) {
  dt <- as.data.table(copy(dt))
  
  for (m in c("Sens","Spec","PPV","NPV")) {
    for (method in c("G1","G12","EvoR")) {
      dt[, paste0(m,"_",method,"_fmt") := fmt_pct_ci(
        safe_col(dt, paste0(m,"_",method)),
        safe_col(dt, paste0(m,"_",method,"_CI_low")),
        safe_col(dt, paste0(m,"_",method,"_CI_high")), d)]
    }
    for (base in c("G1","G12")) {
      dt[, paste0(m,"_Delta_",base,"_fmt") := fmt_pct_ci(
        safe_col(dt, paste0(m,"_Delta_",base)),
        safe_col(dt, paste0(m,"_Delta_",base,"_CI_Boot_low")),
        safe_col(dt, paste0(m,"_Delta_",base,"_CI_Boot_high")), d)]
    }
  }
  
  for (p in c("Sens_P_G1","Sens_P_G12","Spec_P_G1","Spec_P_G12",
              "Sens_P_G1_BH","Sens_P_G12_BH","Spec_P_G1_BH","Spec_P_G12_BH")) {
    dt[, paste0(p,"_fmt") := fmt_p(safe_col(dt, p), p_d)]
  }
  
  # prepend group_col to id columns if provided
  id_cols <- c("Drug", if (!is.null(group_col) && group_col %in% names(dt)) group_col,
               "sample_size", "N_R")
  
  out_cols <- c(
    id_cols,
    paste0(rep(c("Sens","Spec","PPV","NPV"), each=3),
           "_", c("G1","G12","EvoR"), "_fmt"),
    paste0(rep(c("Sens","Spec"), each=4),
           "_Delta_", rep(c("G1","G12"), each=2),
           rep(c("_fmt","_P_fmt","_P_BH_fmt"), times=4)[seq(8)]),
    "PPV_Delta_G1_fmt","PPV_Delta_G12_fmt",
    "NPV_Delta_G1_fmt","NPV_Delta_G12_fmt"
  )
  
  out_cols <- out_cols[out_cols %in% names(dt)]
  dt[, ..out_cols]
}

# ── Read files and run per-drug evaluation ────────────────────────────────────

path_config <- data.table(
  Drug     = drug_list,
  G1_path  = file.path(results_dir, drug_list, "Final_Evaluate", "G1",        "isolate_predictions.tsv"),
  G12_path = file.path(results_dir, drug_list, "Final_Evaluate", "G1_2",      "isolate_predictions.tsv"),
  EvoR_path= file.path(results_dir, drug_list, "Final_Evaluate", "EvoResist", "isolate_predictions.tsv")
)

raw_list       <- list()
raw_strat_list <- list()   # only populated when group_col is set

for (i in seq_len(nrow(path_config))) {
  drug_name <- path_config$Drug[i]
  message("Processing: ", drug_name)
  
  whoG1_dt <- fread(path_config$G1_path[i])
  G12_dt   <- fread(path_config$G12_path[i])
  EvoR_dt  <- fread(path_config$EvoR_path[i])
  
  dt <- merge_prediction_tables(
    whoG1     = whoG1_dt,
    G12       = G12_dt,
    EvoResist = EvoR_dt,
    id_col    = "ID",
    pheno_col = "Phenotype",
    pred_col  = "Prediction",
    drug_col  = NULL,
    group_col = group_col
  )
  dt[, Drug := drug_name]
  
  # Overall (unstratified)
  raw_list[[drug_name]] <- build_result_row(dt, drug_name = drug_name, B = 5000, seed = 20260220)
  
  # Stratified (only when group_col is supplied)
  if (!is.null(group_col)) {
    if (!group_col %in% names(dt)) {
      warning("group_col '", group_col, "' not found in merged table for ",
              drug_name, "; skipping stratification.")
    } else {
      strata <- sort(unique(na.omit(dt[[group_col]])))
      for (s in strata) {
        dt_s <- dt[get(group_col) == s]
        if (nrow(dt_s) == 0) next
        row_s <- build_result_row(dt_s, drug_name = drug_name, B = 5000, seed = 20260220)
        row_s[, (group_col) := s]
        raw_strat_list[[paste(drug_name, s, sep = "_")]] <- row_s
      }
    }
  }
}

# ── Overall: BH-correct, format, write ─────────────────────────────────────────
raw_all <- rbindlist(raw_list, use.names = TRUE, fill = TRUE)
raw_all[, sample_size := N_total]

for (p in c("Sens_P_G1","Sens_P_G12","Spec_P_G1","Spec_P_G12")) {
  raw_all[, paste0(p,"_BH") := p.adjust(get(p), method = "BH")]
}

fwrite(format_result_table(raw_all, group_col = NULL),
       file.path(output_dir, "comparison_formatted.tsv"),
       sep = "	", quote = FALSE, na = "NA")
message("Overall formatted results written.")

# ── Stratified: BH-correct within each stratum, format, write ──────────────────
if (!is.null(group_col) && length(raw_strat_list) > 0) {
  raw_strat <- rbindlist(raw_strat_list, use.names = TRUE, fill = TRUE)
  raw_strat[, sample_size := N_total]
  
  # BH correction within each stratum level across drugs
  for (p in c("Sens_P_G1","Sens_P_G12","Spec_P_G1","Spec_P_G12")) {
    raw_strat[, paste0(p,"_BH") := p.adjust(get(p), method = "BH"), by = group_col]
  }
  
  fwrite(format_result_table(raw_strat, group_col = group_col),
         file.path(output_dir, "comparison_stratified_formatted.tsv"),
         sep = "	", quote = FALSE, na = "NA")
  message("Stratified formatted results written (group_col = '", group_col, "').")
}

message("Done. Output written to ", output_dir)