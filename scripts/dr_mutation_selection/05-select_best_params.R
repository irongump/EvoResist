# 05-select_best_params.R
#
# Two-step parameter selection for the sensitivity analysis pipeline.
#
# --step threshold
#   Reads list2 train/test metrics from the convergence-threshold sweep
#   (promoter = 500 bp, threshold 2-6). Selects the best threshold per
#   drug by maximising MCC on the test set (train MCC as tiebreaker).
#   Writes:
#     <results_dir>/best_thresholds.tsv   – human-readable summary
#     <results_dir>/best_thresholds.sh    – bash-sourceable variable file
#                                           (BEST_THRESHOLD_<DRUG>=<value>)
#                                           sourced by pipeline.sh before
#                                           the promoter sweep
#
# --step promoter
#   Reads list2 train/test metrics from the promoter-length sweep
#   (best threshold per drug, promoter 100-1000 bp). Selects the best
#   promoter length per drug by the same MCC criterion.
#   Writes:
#     <results_dir>/best_promoters.tsv    – human-readable summary
#     <results_dir>/best_promoters.sh     – bash-sourceable variable file
#                                           (BEST_PROMOTER_<DRUG>=<value>)
#
# Usage (called from pipeline.sh):
#   Rscript ./05-select_best_params.R --step threshold --results_dir <path>
#   Rscript ./05-select_best_params.R --step promoter  --results_dir <path>
#
# Optional args:
#   --promoter_fixed <int>   promoter length used during threshold sweep
#                            (default: 500)
#   --metric <string>        metric used for selection: MCC or
#                            Balanced_accuracy (default: MCC)

library(data.table)

# ── Parse arguments ───────────────────────────────────────────────────────────
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out  <- list(
    step           = NULL,
    results_dir    = NULL,
    promoter_fixed = 500L,
    metric         = "MCC"
  )
  i <- 1
  while (i <= length(args)) {
    switch(args[i],
      "--step"           = { out$step           <- args[i + 1]; i <- i + 2 },
      "--results_dir"    = { out$results_dir    <- args[i + 1]; i <- i + 2 },
      "--promoter_fixed" = { out$promoter_fixed <- as.integer(args[i + 1]); i <- i + 2 },
      "--metric"         = { out$metric         <- args[i + 1]; i <- i + 2 },
      { stop("Unknown argument: ", args[i]) }
    )
  }
  if (is.null(out$step))        stop("--step is required (threshold or promoter)")
  if (is.null(out$results_dir)) stop("--results_dir is required")
  if (!out$step %in% c("threshold", "promoter"))
    stop("--step must be 'threshold' or 'promoter'")
  out
}

args <- parse_args()

drug_list <- c("RIF","INH","EMB","PZA","LFX","MFX","BDQ","LZD","AMK","STM","ETO","KAN","CAP")

# ── Helper: compute MCC and balanced accuracy from confusion matrix cols ──────
add_metrics <- function(dt) {
  cols <- c("TP","TN","FP","FN")
  dt[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
  dt[, Balanced_accuracy := (Sensitivity + Specificity) / 2]
  dt[, MCC := (TP * TN - FP * FN) /
               sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))]
  dt
}

# ── Helper: select best parameter combination per drug ────────────────────────
# Selection criterion:
#   1. Maximise test-set metric value
#   2. Tiebreak by maximising train-set metric value
# param_col: name of the column being swept (Threshold or Promoter_length)
select_best <- function(dt, param_col, metric_col) {
  wide <- dcast(
    dt[Set %in% c("Train", "Test")],
    Drug + get(param_col) ~ Set,
    value.var = metric_col
  )
  setnames(wide, "get", param_col)

  best <- wide[
    !is.na(Test),
    .SD[order(-Test, -Train)][1],
    by = Drug
  ]
  setnames(best, c("Test", "Train"),
           c(paste0("best_test_", metric_col), paste0("best_train_", metric_col)))
  setnames(best, param_col, paste0("best_", param_col))
  best
}

# ── Helper: compute overall (train + test combined) sens/spec for best combo ──
overall_metrics <- function(dt, best_dt, param_col) {
  on_cols <- c("Drug", param_col)
  names(on_cols) <- c("Drug", paste0("best_", param_col))

  best_rows <- dt[
    best_dt,
    on = .(Drug),
    allow.cartesian = FALSE,
    nomatch = 0
  ]
  # filter to only the best parameter value
  best_rows <- best_rows[
    mapply(function(d, p) {
      best_val <- best_dt[Drug == d, get(paste0("best_", param_col))]
      p == best_val
    }, Drug, get(param_col))
  ]

  ov <- best_rows[
    , .(TP = sum(TP), TN = sum(TN), FP = sum(FP), FN = sum(FN)),
    by = Drug
  ]
  ov[, `:=`(
    Sensitivity_overall = TP / (TP + FN),
    Specificity_overall = TN / (TN + FP),
    N_R        = TP + FN,
    N_S        = TN + FP,
    N_isolates = TP + TN + FP + FN
  )]
  ov
}

# ── Helper: write bash-sourceable variable file ───────────────────────────────
write_bash_vars <- function(dt, var_prefix, param_col, out_path) {
  lines <- dt[, paste0(var_prefix, "_", Drug, "=", get(paste0("best_", param_col)))]
  writeLines(lines, out_path)
  message("Bash variable file written to ", out_path)
}

# =============================================================================
# STEP: threshold
# Sweep: promoter fixed at args$promoter_fixed, threshold 2-6
# =============================================================================

if (args$step == "threshold") {

  message("Step: threshold selection (promoter = ", args$promoter_fixed, " bp)")

  out_list <- data.table()

  for (drug_name in drug_list) {
    for (threshold in 2:6) {

      train_path <- file.path(
        args$results_dir, drug_name,
        paste0("Threshold_", threshold, "_Promoter_", args$promoter_fixed),
        "train_list2", "overall_metrics.tsv"
      )
      test_path <- file.path(
        args$results_dir, drug_name,
        paste0("Threshold_", threshold, "_Promoter_", args$promoter_fixed),
        "test_list2", "overall_metrics.tsv"
      )

      if (!file.exists(train_path) || !file.exists(test_path)) {
        warning("Missing metrics file for ", drug_name, " threshold=", threshold, "; skipping.")
        next
      }

      train_dt <- fread(train_path)
      train_dt[, `:=`(Drug = drug_name, Threshold = threshold,
                       Promoter_length = args$promoter_fixed, Set = "Train")]

      test_dt <- fread(test_path)
      test_dt[, `:=`(Drug = drug_name, Threshold = threshold,
                      Promoter_length = args$promoter_fixed, Set = "Test")]

      out_list <- rbind(out_list, train_dt, test_dt, fill = TRUE)
    }
  }

  out_list <- add_metrics(out_list)

  best <- select_best(out_list, param_col = "Threshold", metric_col = args$metric)
  ov   <- overall_metrics(out_list, best, param_col = "Threshold")
  result <- merge(best, ov, by = "Drug")

  # Write TSV summary
  tsv_path <- file.path(args$results_dir, "best_thresholds.tsv")
  fwrite(result, tsv_path, sep = "\t", col.names = TRUE)
  message("Best thresholds written to ", tsv_path)

  # Write bash-sourceable file
  sh_path <- file.path(args$results_dir, "best_thresholds.sh")
  write_bash_vars(result, var_prefix = "BEST_THRESHOLD", param_col = "Threshold", out_path = sh_path)

  # Print summary to stdout
  print(result[order(Drug), .(Drug, best_Threshold,
                               Sensitivity_overall, Specificity_overall)])
}

# =============================================================================
# STEP: promoter
# Sweep: best threshold per drug (from best_thresholds.tsv), promoter 100-1000
# =============================================================================

if (args$step == "promoter") {

  message("Step: promoter length selection")

  # Load best thresholds determined in the previous step
  thres_path <- file.path(args$results_dir, "best_thresholds.tsv")
  if (!file.exists(thres_path))
    stop("best_thresholds.tsv not found at ", thres_path,
         ". Run --step threshold first.")

  best_thresholds <- fread(thres_path)[, .(Drug, best_Threshold)]

  out_list <- data.table()

  for (drug_name in drug_list) {
    threshold <- best_thresholds[Drug == drug_name, best_Threshold]

    for (promoter_region in seq(100, 1000, by = 100)) {

      train_path <- file.path(
        args$results_dir, drug_name, "Sens_analysis_75",
        paste0("Threshold_", threshold, "_Promoter_", promoter_region),
        "train_list2", "overall_metrics.tsv"
      )
      test_path <- file.path(
        args$results_dir, drug_name, "Sens_analysis_75",
        paste0("Threshold_", threshold, "_Promoter_", promoter_region),
        "test_list2", "overall_metrics.tsv"
      )

      if (!file.exists(train_path) || !file.exists(test_path)) {
        warning("Missing metrics file for ", drug_name,
                " threshold=", threshold, " promoter=", promoter_region, "; skipping.")
        next
      }

      train_dt <- fread(train_path)
      train_dt[, `:=`(Drug = drug_name, Threshold = threshold,
                       Promoter_length = promoter_region, Set = "Train")]

      test_dt <- fread(test_path)
      test_dt[, `:=`(Drug = drug_name, Threshold = threshold,
                      Promoter_length = promoter_region, Set = "Test")]

      out_list <- rbind(out_list, train_dt, test_dt, fill = TRUE)
    }
  }

  out_list <- add_metrics(out_list)

  best <- select_best(out_list, param_col = "Promoter_length", metric_col = args$metric)
  ov   <- overall_metrics(out_list, best, param_col = "Promoter_length")
  result <- merge(best, ov, by = "Drug")
  result <- merge(result, best_thresholds, by = "Drug")

  # Write TSV summary
  tsv_path <- file.path(args$results_dir, "best_promoters.tsv")
  fwrite(result, tsv_path, sep = "\t", col.names = TRUE)
  message("Best promoters written to ", tsv_path)

  # Write bash-sourceable file
  sh_path <- file.path(args$results_dir, "best_promoters.sh")
  write_bash_vars(result, var_prefix = "BEST_PROMOTER",
                  param_col = "Promoter_length", out_path = sh_path)

  # Print summary to stdout
  print(result[order(Drug), .(Drug, best_Threshold, best_Promoter_length,
                               Sensitivity_overall, Specificity_overall)])
}