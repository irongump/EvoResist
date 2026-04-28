#!/usr/bin/env python3
"""
06-gain_evaluation.py

Incremental gain ranking of EvoResist candidates relative to the shared part with WHO G1+G2.

Baseline definition:
    intersection of WHO G1+G2 and EvoResist,
    i.e. variants present in both lists.

Candidates:
    variants in EvoResist that are absent from WHO G1+G2.

For each candidate, computes the marginal improvement in sensitivity,
specificity, and PPV when that variant is added individually to the baseline.
The "full panel" reference for cumulative gain is List2 alone.

Variants are ranked by (Delta_Sens, Delta_Spec, IncR, Delta_PPV) descending,
and cumulative sensitivity gain is tracked across the ranked list.

Output: a single TSV file with per-variant rankings and cumulative metrics.
"""

import pandas as pd
import argparse
import os
import sys
import math
from collections import defaultdict

# ==========================================================
# CLI
# ==========================================================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Incremental gain ranking using List1∩List2 as baseline and List2-List1 as candidates."
    )

    parser.add_argument('--list1_variants_file', required=True,
                        help="Baseline variant list (e.g. WHO G1+G2).")
    parser.add_argument('--list2_variants_file', required=True,
                        help="Candidate variant list (e.g. EvoResist final list).")
    parser.add_argument('--id_list_file', required=True,
                        help="Per-drug sample list TSV with columns Run and pheno.")

    parser.add_argument('--snp_dir', required=True,
                        help="Directory containing per-sample SNP annotation files ({ID}.ano).")
    parser.add_argument('--additional_snp_dir', default=None,
                        help="Optional fallback SNP directory.")

    parser.add_argument('--indel_dir', required=True,
                        help="Directory containing per-sample indel annotation files ({ID}.indel.ano).")
    parser.add_argument('--additional_indel_dir', default=None,
                        help="Optional fallback indel directory.")

    parser.add_argument('--output_file', default='incremental_gain_ranking.tsv',
                        help="Output TSV file path.")

    return parser.parse_args()

# ==========================================================
# Wilson CI (kept for consistency, not used in ranking)
# ==========================================================

def wilson_ci(x, n, z=1.96):
    if n == 0:
        return (None, None)
    phat = x / n
    denom = 1 + z**2 / n
    center = (phat + z**2/(2*n)) / denom
    half = (z * math.sqrt((phat*(1-phat)/n) + z**2/(4*n*n))) / denom
    return (max(0, center-half), min(1, center+half))

# ==========================================================
# Metric computation
# ==========================================================

def compute_metrics(pred_R, actual_R, actual_S, universe):
    """Return (sensitivity, specificity, PPV) for given predictions."""
    pred_R = pred_R & universe
    pred_S = universe - pred_R

    TP = len(pred_R & actual_R)
    FP = len(pred_R & actual_S)
    TN = len(pred_S & actual_S)
    FN = len(pred_S & actual_R)

    sens = TP/(TP+FN) if TP+FN>0 else 0
    spec = TN/(TN+FP) if TN+FP>0 else 0
    ppv  = TP/(TP+FP) if TP+FP>0 else 0

    return sens, spec, ppv

# ==========================================================
# Load ID list
# ==========================================================

def load_id_list_new(file):
    """
    Load per-drug sample list. Expects columns Run and pheno.
    IDs with inconsistent phenotypes across duplicate rows are excluded.
    """
    id_df = pd.read_csv(file, sep="\t", dtype=str)
    id_df = id_df.rename(columns={"Run":"ID", "pheno":"Phenotype"})

    id_df.dropna(subset=["ID","Phenotype"], inplace=True)
    id_df["ID"]        = id_df["ID"].str.strip()
    id_df["Phenotype"] = id_df["Phenotype"].str.strip().str.upper()

    duplicated   = id_df[id_df.duplicated(subset=["ID"], keep=False)]
    pheno_counts = duplicated.groupby("ID")["Phenotype"].nunique()
    inconsistent = pheno_counts[pheno_counts>1].index.tolist()

    if inconsistent:
        id_df = id_df[~id_df["ID"].isin(inconsistent)]
        print(f"Excluded {len(inconsistent)} inconsistent IDs")

    id_df = id_df.drop_duplicates(subset=["ID"], keep="first")
    return id_df

# ==========================================================
# Load SNP / indel files (dual directory support)
# ==========================================================

def load_all_snp_files(id_df, snp_dir, additional_snp_dir=None):
    """
    Load SNP annotation files ({ID}.ano) for all IDs.
    Falls back to additional_snp_dir if the primary file is absent.
    Returns a dict mapping ID -> set of (position, ref, alt) tuples.
    """
    snp_dict = {}

    for _, row in id_df.iterrows():
        id_       = row["ID"]
        primary   = os.path.join(snp_dir, f"{id_}.ano")
        secondary = os.path.join(additional_snp_dir, f"{id_}.ano") if additional_snp_dir else None

        if os.path.isfile(primary):
            snp_file = primary
        elif secondary and os.path.isfile(secondary):
            snp_file = secondary
        else:
            continue

        snp_df = pd.read_csv(snp_file, sep="\t", header=None,
                             names=["position","ref","alt"],
                             usecols=[0,1,2], dtype=str)
        snp_df["position"] = snp_df["position"].astype(str).str.strip()
        snp_df["ref"]      = snp_df["ref"].astype(str).str.upper().str.strip()
        snp_df["alt"]      = snp_df["alt"].astype(str).str.upper().str.strip()
        snp_dict[id_]      = set(zip(snp_df["position"], snp_df["ref"], snp_df["alt"]))

    return snp_dict


def load_all_indel_files(id_df, indel_dir, additional_indel_dir=None):
    """
    Load indel annotation files ({ID}.indel.ano) for all IDs.
    Falls back to additional_indel_dir if the primary file is absent.
    Returns a dict mapping ID -> set of (position, ref, alt) tuples.
    """
    indel_dict = {}

    for _, row in id_df.iterrows():
        id_       = row["ID"]
        primary   = os.path.join(indel_dir, f"{id_}.indel.ano")
        secondary = os.path.join(additional_indel_dir, f"{id_}.indel.ann") if additional_indel_dir else None

        if os.path.isfile(primary):
            indel_file = primary
        elif secondary and os.path.isfile(secondary):
            indel_file = secondary
        else:
            continue

        indel_df = pd.read_csv(indel_file, sep="\t", header=None,
                               names=["position","ref","alt"],
                               usecols=[0,1,2], dtype=str)
        indel_df["position"] = indel_df["position"].astype(str).str.strip()
        indel_df["ref"]      = indel_df["ref"].astype(str).str.upper().str.strip()
        indel_df["alt"]      = indel_df["alt"].astype(str).str.upper().str.strip()
        indel_dict[id_]      = set(zip(indel_df["position"], indel_df["ref"], indel_df["alt"]))

    return indel_dict

# ==========================================================
# Variant list loader (supports headered and legacy headerless formats)
# ==========================================================

def load_variant_list(filepath, label):
    """
    Load a variant list TSV.

    Auto-detects format:
      - Headered (new): first row contains 'position', 'ref', 'alt' as column names.
      - Headerless (legacy): fixed column order [variant, position, ref, alt, grading].

    Only position, ref, and alt are used downstream.
    """
    with open(filepath) as fh:
        first_line_cols = [c.strip() for c in fh.readline().rstrip("\n").split("\t")]

    required = {"position", "ref", "alt"}

    if required.issubset(set(first_line_cols)):
        df = pd.read_csv(filepath, sep="\t", dtype=str)
        df.columns = [c.strip() for c in df.columns]
    else:
        df = pd.read_csv(filepath, sep="\t", header=None,
                         names=["variant", "position", "ref", "alt", "grading"],
                         dtype=str)

    missing_cols = required - set(df.columns)
    if missing_cols:
        print(f"ERROR: {label} file is missing required columns: {sorted(missing_cols)}", file=sys.stderr)
        sys.exit(1)

    df.dropna(subset=["position", "ref", "alt"], inplace=True)
    df["position"] = df["position"].astype(str).str.strip()
    df["ref"]      = df["ref"].astype(str).str.upper().str.strip()
    df["alt"]      = df["alt"].astype(str).str.upper().str.strip()

    print(f"Loaded {len(df)} variants from {label}.")
    return df

# ==========================================================
# Main
# ==========================================================

def main():

    args = parse_arguments()

    # Load variant lists
    list1_df  = load_variant_list(args.list1_variants_file, "List1")
    list2_df  = load_variant_list(args.list2_variants_file, "List2")
    list1_set = set(zip(list1_df["position"], list1_df["ref"], list1_df["alt"]))
    list2_set = set(zip(list2_df["position"], list2_df["ref"], list2_df["alt"]))

    # Load IDs and genotype files
    id_df      = load_id_list_new(args.id_list_file)
    snp_dict   = load_all_snp_files(id_df, args.snp_dir, args.additional_snp_dir)
    indel_dict = load_all_indel_files(id_df, args.indel_dir, args.additional_indel_dir)

    included_ids = set(snp_dict.keys()) & set(indel_dict.keys())
    id_df        = id_df[id_df["ID"].isin(included_ids)]

    universe = set(id_df["ID"])
    actual_R = set(id_df[id_df["Phenotype"] == "R"]["ID"])
    actual_S = set(id_df[id_df["Phenotype"] == "S"]["ID"])

    # Build variant -> ID presence mapping
    variant_to_ids = defaultdict(set)
    for id_ in universe:
        for v in snp_dict[id_] | indel_dict[id_]:
            variant_to_ids[v].add(id_)

    # ------------------------------------------------------------------
    # Baseline: List1 ∩ List2 (variants shared by both lists).
    # This isolates the contribution of EvoResist-only variants from
    # the background of variants that both catalogues agree on.
    # ------------------------------------------------------------------
    baseline_set  = list1_set & list2_set
    baseline_pred = set()
    for v in baseline_set:
        baseline_pred |= variant_to_ids.get(v, set())

    S0, Spec0, PPV0 = compute_metrics(baseline_pred, actual_R, actual_S, universe)

    # ------------------------------------------------------------------
    # Full panel: List2 only.
    # Used as the upper bound for cumulative sensitivity gain.
    # ------------------------------------------------------------------
    full_pred = set()
    for v in list2_set:
        full_pred |= variant_to_ids.get(v, set())

    S_all, _, _ = compute_metrics(full_pred, actual_R, actual_S, universe)
    total_gain   = S_all - S0

    # ------------------------------------------------------------------
    # Candidates: List2 - List1.
    # For each candidate, compute marginal metrics when added to baseline.
    # ------------------------------------------------------------------
    candidates = list2_set - list1_set
    rows       = []

    for v in candidates:
        pos, ref, alt = v
        new_pred      = baseline_pred | variant_to_ids.get(v, set())
        S_v, Spec_v, PPV_v = compute_metrics(new_pred, actual_R, actual_S, universe)

        # Number of R isolates newly called by this variant beyond the baseline
        inc_R = len((variant_to_ids.get(v, set()) & actual_R) - (baseline_pred & actual_R))

        rows.append({
            "Position":   pos,
            "Ref":        ref,
            "Alt":        alt,
            "Delta_Sens": S_v    - S0,
            "Delta_Spec": Spec_v - Spec0,
            "Delta_PPV":  PPV_v  - PPV0,
            "IncR":       inc_R,
        })

    df = pd.DataFrame(rows)

    # Rank by (Delta_Sens, Delta_Spec, IncR, Delta_PPV) descending
    df = df.sort_values(
        by=["Delta_Sens", "Delta_Spec", "IncR", "Delta_PPV"],
        ascending=[False, False, False, False]
    ).reset_index(drop=True)
    df["Rank"] = df.index + 1

    # ------------------------------------------------------------------
    # Cumulative sensitivity gain in ranked order
    # ------------------------------------------------------------------
    cum_pred = baseline_pred.copy()
    sens_cum = []
    gain_cum = []
    prop_cum = []

    for _, row in df.iterrows():
        v = (row["Position"], row["Ref"], row["Alt"])
        cum_pred |= variant_to_ids.get(v, set())
        S_cum, _, _ = compute_metrics(cum_pred, actual_R, actual_S, universe)
        G_cum = S_cum - S0
        sens_cum.append(S_cum)
        gain_cum.append(G_cum)
        prop_cum.append(G_cum / total_gain if total_gain > 0 else 0)

    df["Sens_cum"] = sens_cum
    df["Cum_Gain"] = gain_cum
    df["Cum_Prop"] = prop_cum

    df.to_csv(args.output_file, sep="\t", index=False)
    print("Incremental gain ranking saved to:", args.output_file)


if __name__ == "__main__":
    main()
