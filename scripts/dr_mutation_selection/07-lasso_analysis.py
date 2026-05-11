#!/usr/bin/env python3
"""
07-lasso_analysis.py

Penalized logistic regression (L1 or L2) for variant-phenotype association,
with lineage covariates and variant × lineage interaction terms.

Model structure:
    X = [variant presence (0/1)] | [lineage dummies, ref=LINEAGE4]
                                  | [variant × lineage interactions]

Variants present in fewer than MIN_ISOLATES isolates are dropped before
model fitting. Interaction columns are derived from the filtered variant
set only (strong-hierarchy principle).

After L1 (LASSO) fitting, an unpenalized refit is performed on the
selected (non-zero) features using statsmodels.Logit, yielding standard
errors, p-values, and 95% CIs. Falls back to sklearn penalty=None when
statsmodels is not installed.

Fixed hyperparameters (edit at top of script if needed):
    MIN_ISOLATES = 50          cohort-wide MAC filter threshold
    CV_FOLDS     = 5           cross-validation folds
    MAX_ITER     = 5000        solver iteration limit
    CS           = logspace(-4, 4, 15)   regularization grid

Usage:
    python3 07-lasso_analysis.py \\
        --variants_file <path> \\
        --id_list_file  <path> \\
        --snp_dir       <path> \\
        --indel_dir     <path> \\
        --penalty       l1|l2  \\
        --output_dir    <path>

Output files (written to --output_dir):
    data_summary.tsv          cohort and feature-matrix summary
    cv_performance.tsv        cross-validated AUC across C grid
    model_coefficients.tsv    penalized model coefficients and ORs
    refit_coefficients.tsv    unpenalized refit on LASSO-selected features
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd

from scipy.sparse import csr_matrix, hstack
from sklearn.linear_model import LogisticRegressionCV
from sklearn.exceptions import ConvergenceWarning
import warnings

try:
    import statsmodels.api as sm
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False

# ==========================================================
# Fixed hyperparameters
# ==========================================================

MIN_ISOLATES = 50
CV_FOLDS     = 5
MAX_ITER     = 5000
CS           = np.logspace(-4, 4, 15)

# ==========================================================
# Constants
# ==========================================================

VALID_LINEAGES = [f"LINEAGE{i}" for i in range(1, 8)]
REF_LINEAGE    = "LINEAGE4"

# ==========================================================
# CLI
# ==========================================================

def parse_arguments():
    p = argparse.ArgumentParser(
        description="Penalized logistic regression with variants, lineage dummies, "
                    "and variant×lineage interaction terms."
    )
    p.add_argument('--variants_file', required=True,
                   help="Variant list TSV with header. Required columns: position, ref, alt.")
    p.add_argument('--id_list_file', required=True,
                   help="Sample list TSV. Required columns: Run, pheno, Lineage.")
    p.add_argument('--snp_dir', required=True,
                   help="Directory containing SNP annotation files ({ID}.ano).")
    p.add_argument('--indel_dir', required=True,
                   help="Directory containing indel annotation files ({ID}.indel.ano).")
    p.add_argument('--penalty', required=True, choices=['l1', 'l2'],
                   help="Regularization penalty: l1 (LASSO) or l2 (ridge).")
    p.add_argument('--output_dir', required=True,
                   help="Output directory.")
    return p.parse_args()

# ==========================================================
# Data loading
# ==========================================================

def load_variant_list(variants_file):
    """
    Load variant list TSV. Auto-detects whether a header is present.
    New format: header row containing at least 'position', 'ref', 'alt'.
    Legacy format: no header, fixed column order [variant, position, ref, alt, grading].
    Only position, ref, and alt are used downstream.
    """
    with open(variants_file) as fh:
        first_cols = [c.strip() for c in fh.readline().rstrip("\n").split("\t")]

    required = {"position", "ref", "alt"}
    if required.issubset(set(first_cols)):
        df = pd.read_csv(variants_file, sep='\t', dtype=str, engine='python')
        df.columns = [c.strip() for c in df.columns]
    else:
        df = pd.read_csv(variants_file, sep='\t', header=None,
                         names=['variant', 'position', 'ref', 'alt', 'grading'],
                         dtype=str, engine='python')

    missing_cols = required - set(df.columns)
    if missing_cols:
        raise ValueError(f"Variant file missing required columns: {sorted(missing_cols)}")

    df.dropna(subset=['position', 'ref', 'alt'], inplace=True)
    df['position'] = df['position'].astype(str).str.strip()
    df['ref']      = df['ref'].astype(str).str.upper().str.strip()
    df['alt']      = df['alt'].astype(str).str.upper().str.strip()
    return df


def load_id_list(id_list_file):
    """
    Load per-drug sample list. Required columns: Run, pheno, Lineage.
    Samples are restricted to LINEAGE1–LINEAGE7; all others are excluded.
    IDs with conflicting phenotype entries are removed before returning.
    """
    df = pd.read_csv(id_list_file, sep='\t', dtype=str, low_memory=False)
    for col in ('Run', 'pheno', 'Lineage'):
        if col not in df.columns:
            raise ValueError(f"ID file must contain column '{col}'.")

    df = df.rename(columns={'Run': 'ID', 'pheno': 'Phenotype'})
    df = df[['ID', 'Phenotype', 'Lineage']].copy()
    df.dropna(subset=['ID', 'Phenotype', 'Lineage'], inplace=True)
    df['ID']        = df['ID'].astype(str).str.strip()
    df['Phenotype'] = df['Phenotype'].astype(str).str.strip().str.upper()
    df['Lineage']   = df['Lineage'].astype(str).str.strip().str.upper()

    bad_pheno = set(df['Phenotype']) - {'R', 'S'}
    if bad_pheno:
        raise ValueError(f"Invalid phenotypes: {bad_pheno}. Only 'R' or 'S' are allowed.")

    df = df[df['Lineage'].isin(VALID_LINEAGES)].copy()

    dups = df[df.duplicated(subset=['ID'], keep=False)]
    if not dups.empty:
        inconsistent = dups.groupby('ID')['Phenotype'].nunique()
        inconsistent = inconsistent[inconsistent > 1].index.tolist()
        if inconsistent:
            df = df[~df['ID'].isin(inconsistent)].copy()
        df = df.drop_duplicates(subset=['ID'], keep='first')

    return df


def load_snp_files(id_df, snp_dir):
    """
    Load SNP annotation files ({ID}.ano) for all IDs.
    Returns a dict mapping ID -> set of (position, ref, alt) tuples.
    """
    snp_dict = {}
    missing  = 0
    total    = len(id_df)

    for i, (_, row) in enumerate(id_df.iterrows()):
        id_  = row["ID"]
        path = os.path.join(snp_dir, f"{id_}.ano")
        if not os.path.isfile(path):
            missing += 1
            continue
        try:
            df = pd.read_csv(path, sep="\t", header=None,
                             names=["position", "ref", "alt"],
                             usecols=[0, 1, 2], dtype=str)
            df.dropna(subset=["position", "ref", "alt"], inplace=True)
            if df.empty:
                missing += 1
                continue
            df["position"] = df["position"].astype(str).str.strip()
            df["ref"]      = df["ref"].str.upper().str.strip()
            df["alt"]      = df["alt"].str.upper().str.strip()
            snp_dict[id_]  = set(zip(df["position"], df["ref"], df["alt"]))
        except Exception as e:
            print(f"Warning: could not read {path}: {e}")
            missing += 1

        if (i + 1) % 1000 == 0 or (i + 1) == total:
            print(f"  SNP: processed {i+1}/{total}")

    print(f"SNP files loaded: {len(snp_dict)} IDs. Missing: {missing}.")
    return snp_dict, missing


def load_indel_files(id_df, indel_dir):
    """
    Load indel annotation files ({ID}.indel.ano) for all IDs.
    Returns a dict mapping ID -> set of (position, ref, alt) tuples.
    """
    indel_dict = {}
    missing    = 0
    total      = len(id_df)

    for i, (_, row) in enumerate(id_df.iterrows()):
        id_  = row["ID"]
        path = os.path.join(indel_dir, f"{id_}.indel.ano")
        if not os.path.isfile(path):
            missing += 1
            continue
        try:
            df = pd.read_csv(path, sep="\t", header=None,
                             names=["position", "ref", "alt"],
                             usecols=[0, 1, 2], dtype=str)
            df.dropna(subset=["position", "ref", "alt"], inplace=True)
            if df.empty:
                missing += 1
                continue
            df["position"]   = df["position"].astype(str).str.strip()
            df["ref"]        = df["ref"].str.upper().str.strip()
            df["alt"]        = df["alt"].str.upper().str.strip()
            indel_dict[id_]  = set(zip(df["position"], df["ref"], df["alt"]))
        except Exception as e:
            print(f"Warning: could not read {path}: {e}")
            missing += 1

        if (i + 1) % 1000 == 0 or (i + 1) == total:
            print(f"  INDEL: processed {i+1}/{total}")

    print(f"INDEL files loaded: {len(indel_dict)} IDs. Missing: {missing}.")
    return indel_dict, missing

# ==========================================================
# Design matrix
# ==========================================================

def build_design_matrix(id_df, snp_dict, indel_dict, variant_keys):
    """
    Build a sparse design matrix:
        X = [variants (0/1)] | [lineage dummies, ref=LINEAGE4]
                             | [variant × lineage interactions]

    Variant columns are filtered to those present in >= MIN_ISOLATES isolates.
    Interaction columns are built from the filtered variant set only.

    Returns:
        X                  csr_matrix, shape (n_samples, n_features)
        y                  np.ndarray of 0/1 (1=R)
        feature_names      list of str
        ids                list of sample IDs in row order
        kept_variant_keys  list of (pos, ref, alt) after MAC filter
        lineage_levels     list of non-reference lineage labels used
    """
    included_ids = sorted(set(snp_dict.keys()) & set(indel_dict.keys()))
    id_df = id_df[id_df['ID'].isin(included_ids)].set_index('ID').loc[included_ids].reset_index()

    ids = id_df['ID'].tolist()
    n   = len(ids)
    y   = (id_df['Phenotype'].values == 'R').astype(int)

    lineage        = id_df['Lineage'].values
    lineage_levels = [L for L in VALID_LINEAGES if L != REF_LINEAGE]
    lineage_index  = {L: j for j, L in enumerate(lineage_levels)}
    k              = len(lineage_levels)

    combined  = {id_: snp_dict[id_] | indel_dict[id_] for id_ in ids}
    p         = len(variant_keys)
    var_index = {v: j for j, v in enumerate(variant_keys)}

    row_idx, col_idx = [], []
    for i, id_ in enumerate(ids):
        for v in combined[id_]:
            j = var_index.get(v)
            if j is not None:
                row_idx.append(i)
                col_idx.append(j)

    X_var = csr_matrix(
        (np.ones(len(row_idx), dtype=np.float32), (row_idx, col_idx)),
        shape=(n, p), dtype=np.float32
    )

    # Cohort-wide MAC filter
    col_sums          = np.asarray(X_var.sum(axis=0)).ravel()
    keep_cols         = np.where(col_sums >= MIN_ISOLATES)[0]
    X_var             = X_var[:, keep_cols]
    kept_variant_keys = [variant_keys[j] for j in keep_cols]
    print(f"Variant MAC filter: kept {len(kept_variant_keys)}/{p} variants "
          f"with >= {MIN_ISOLATES} isolates.")

    # Lineage dummy matrix (reference category dropped)
    lin_rows, lin_cols = [], []
    for i, L in enumerate(lineage):
        j = lineage_index.get(L)
        if j is not None:
            lin_rows.append(i)
            lin_cols.append(j)
    X_lin = csr_matrix(
        (np.ones(len(lin_rows), dtype=np.float32), (lin_rows, lin_cols)),
        shape=(n, k), dtype=np.float32
    )

    # Variant × lineage interaction blocks
    variant_names      = [f"VAR:{pos}:{ref}:{alt}" for (pos, ref, alt) in kept_variant_keys]
    lineage_names      = [f"LIN:{L}" for L in lineage_levels]
    interaction_blocks = []
    interaction_names  = []

    for j, L in enumerate(lineage_levels):
        block = X_var.multiply(X_lin[:, j])
        interaction_blocks.append(block)
        for vn in variant_names:
            interaction_names.append(f"{vn}*{lineage_names[j]}")

    X_int = hstack(interaction_blocks, format='csr') if interaction_blocks \
            else csr_matrix((n, 0), dtype=np.float32)

    X             = hstack([X_var, X_lin, X_int], format='csr')
    feature_names = variant_names + lineage_names + interaction_names

    return X, y, feature_names, ids, kept_variant_keys, lineage_levels

# ==========================================================
# Model fitting
# ==========================================================

def fit_penalized_logreg(X, y, penalty):
    """
    Fit LogisticRegressionCV (saga solver, balanced class weights, AUC scoring).
    Returns the fitted model and a CV performance table (C vs mean/std AUC).
    """
    warnings.filterwarnings("ignore", category=ConvergenceWarning)

    model = LogisticRegressionCV(
        Cs=CS, cv=CV_FOLDS, penalty=penalty, solver='saga',
        scoring='roc_auc', class_weight='balanced',
        max_iter=MAX_ITER, n_jobs=-1, refit=True
    )
    model.fit(X, y)

    scores   = model.scores_[1] if isinstance(model.scores_, dict) else model.scores_
    cv_df    = pd.DataFrame({
        "C":        model.Cs_,
        "mean_auc": scores.mean(axis=0),
        "std_auc":  scores.std(axis=0),
        "best_C":   float(model.C_[0])
    })
    return model, cv_df


def write_coefficients(output_path, feature_names, coef):
    """Write penalized model coefficients, absolute values, ORs, and nonzero flag."""
    pd.DataFrame({
        "Feature": feature_names,
        "Coef":    coef,
        "AbsCoef": np.abs(coef),
        "OR":      np.exp(coef),
        "NonZero": (coef != 0).astype(int)
    }).sort_values("AbsCoef", ascending=False).to_csv(output_path, sep='\t', index=False)

# ==========================================================
# Unpenalized refit on LASSO-selected features
# ==========================================================

def refit_unpenalized(X, y, feature_names, coef_penalized, output_path):
    """
    Refit an unpenalized logistic regression on the features selected by LASSO
    (non-zero coefficients) to obtain unbiased estimates, SEs, p-values, and CIs.
    Uses statsmodels.Logit; falls back to sklearn penalty=None if not installed.
    SE > 10 is flagged as possible complete separation.
    """
    nonzero_idx = np.where(coef_penalized != 0)[0]
    n_selected  = len(nonzero_idx)

    if n_selected == 0:
        print("Refit skipped: LASSO selected 0 features.")
        pd.DataFrame(columns=[
            "Feature", "Coef", "SE", "z", "pvalue",
            "CI_lower_95", "CI_upper_95", "OR", "OR_CI_lower_95", "OR_CI_upper_95"
        ]).to_csv(output_path, sep='\t', index=False)
        return

    print(f"Refit: {n_selected} features selected → fitting unpenalized logistic regression.")
    X_sel          = X[:, nonzero_idx]
    selected_names = [feature_names[j] for j in nonzero_idx]

    if HAS_STATSMODELS:
        X_dense      = np.asarray(X_sel.todense(), dtype=np.float64)
        X_with_const = sm.add_constant(X_dense, has_constant='add')
        try:
            result = sm.Logit(y.astype(np.float64), X_with_const).fit(
                method='newton', maxiter=MAX_ITER, disp=False, warn_convergence=True)
        except Exception as e:
            print(f"  Newton failed ({e}), retrying with bfgs.")
            try:
                result = sm.Logit(y.astype(np.float64), X_with_const).fit(
                    method='bfgs', maxiter=MAX_ITER, disp=False, warn_convergence=True)
            except Exception as e2:
                print(f"  Refit failed: {e2}. Skipping.")
                return

        params, bse = result.params[1:], result.bse[1:]
        ci          = result.conf_int(alpha=0.05)
        ci_lower, ci_upper = ci[1:, 0], ci[1:, 1]

        out_df = pd.DataFrame({
            "Feature":        selected_names,
            "Coef":           params,
            "SE":             bse,
            "z":              result.tvalues[1:],
            "pvalue":         result.pvalues[1:],
            "CI_lower_95":    ci_lower,
            "CI_upper_95":    ci_upper,
            "OR":             np.exp(params),
            "OR_CI_lower_95": np.exp(ci_lower),
            "OR_CI_upper_95": np.exp(ci_upper),
            "Warning":        np.where(bse > 10, "large_SE_possible_separation", "")
        })
    else:
        print("  statsmodels not found; using sklearn penalty=None (no SE/p-value).")
        from sklearn.linear_model import LogisticRegression
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        lr = LogisticRegression(penalty=None, solver='lbfgs',
                                max_iter=MAX_ITER, fit_intercept=True)
        lr.fit(X_sel, y)
        coef_refit = lr.coef_.ravel()
        out_df = pd.DataFrame({
            "Feature":        selected_names,
            "Coef":           coef_refit,
            "SE":             np.nan, "z": np.nan, "pvalue": np.nan,
            "CI_lower_95":    np.nan, "CI_upper_95": np.nan,
            "OR":             np.exp(coef_refit),
            "OR_CI_lower_95": np.nan, "OR_CI_upper_95": np.nan,
            "Warning":        "statsmodels_not_installed"
        })

    out_df.sort_values("pvalue", na_position='last').reset_index(drop=True)\
          .to_csv(output_path, sep='\t', index=False)
    print(f"  Refit results saved to: {output_path}")

# ==========================================================
# Main
# ==========================================================

def main():
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)

    variants_df  = load_variant_list(args.variants_file)
    variant_keys = sorted(set(zip(variants_df['position'], variants_df['ref'], variants_df['alt'])))
    id_df        = load_id_list(args.id_list_file)

    snp_dict,   missing_snp   = load_snp_files(id_df, args.snp_dir)
    indel_dict, missing_indel = load_indel_files(id_df, args.indel_dir)

    id_df = id_df[id_df['ID'].isin(set(snp_dict.keys()) & set(indel_dict.keys()))].copy()

    X, y, feature_names, ids, kept_variant_keys, lineage_levels = build_design_matrix(
        id_df, snp_dict, indel_dict, variant_keys
    )

    if len(y) == 0:
        print("No isolates remaining after filtering. Exiting.", file=sys.stderr)
        sys.exit(1)

    n, nR, nS = len(y), int(y.sum()), int((y == 0).sum())

    pd.DataFrame([{
        "n_isolates":                   n,
        "n_R":                          nR,
        "n_S":                          nS,
        "n_variants_input":             len(variant_keys),
        "min_isolates_threshold":       MIN_ISOLATES,
        "n_variants_kept_after_filter": len(kept_variant_keys),
        "n_lineage_levels_used":        len(lineage_levels) + 1,
        "ref_lineage":                  REF_LINEAGE,
        "missing_snp_files":            missing_snp,
        "missing_indel_files":          missing_indel,
        "penalty":                      args.penalty,
        "cv_folds":                     CV_FOLDS,
        "X_shape_rows":                 X.shape[0],
        "X_shape_cols":                 X.shape[1]
    }]).to_csv(os.path.join(args.output_dir, "data_summary.tsv"), sep="\t", index=False)

    model, cv_df = fit_penalized_logreg(X, y, penalty=args.penalty)
    cv_df.to_csv(os.path.join(args.output_dir, "cv_performance.tsv"), sep="\t", index=False)

    coef = model.coef_.ravel()
    write_coefficients(os.path.join(args.output_dir, "model_coefficients.tsv"), feature_names, coef)

    print(f"\nDone. Output: {args.output_dir}")
    print(f"Isolates: {n} (R={nR}, S={nS}), Feature matrix: {X.shape}, Best C: {float(model.C_[0])}")
    if args.penalty == 'l1':
        print(f"Nonzero coefficients: {int((coef != 0).sum())}")

    refit_unpenalized(
        X=X, y=y,
        feature_names=feature_names,
        coef_penalized=coef,
        output_path=os.path.join(args.output_dir, "refit_coefficients.tsv"),
    )


if __name__ == "__main__":
    main()
