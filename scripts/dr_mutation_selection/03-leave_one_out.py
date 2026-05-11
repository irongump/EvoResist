#!/usr/bin/env python3
"""
Leave-one-out variant list evaluation.

For a given variant list and ID list (with phenotype), computes:
  1. Overall sensitivity, specificity, precision (PPV), and NPV with 95% Wilson CIs.
  2. Per-variant R/S counts and leave-one-out delta metrics (how each variant
     individually affects each performance metric when removed).
  3. Per-isolate resistance predictions.

Output files written to --output_dir:
  overall_metrics.tsv
  per_variant_analysis.tsv
  isolate_predictions.tsv
"""

import pandas as pd
import argparse
import os
import sys
import math
from collections import defaultdict


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Single variant list: overall metrics with CI + per-variant leave-one-out analysis."
    )
    parser.add_argument("--variants_file",  required=True,
                        help="Variant list TSV with header. Required columns: variant, position, ref, alt.")
    parser.add_argument("--id_list_file",   required=True,
                        help="ID list TSV with columns Run and pheno.")
    parser.add_argument("--snp_dir",        required=True,
                        help="Directory containing SNP annotation files named {ID}.ano.")
    parser.add_argument("--indel_dir",      required=True,
                        help="Directory containing indel annotation files named {ID}.indel.ano.")
    parser.add_argument("--output_dir",     default="variant_analysis_output",
                        help="Directory to write output files.")
    return parser.parse_args()


# ------------------------------------------------------------------
# Statistical helpers
# ------------------------------------------------------------------

def wilson_ci(x, n, z=1.96):
    """Wilson score interval. Returns (lower, upper); (None, None) if n == 0."""
    if n <= 0:
        return (None, None)
    phat   = x / n
    denom  = 1 + z ** 2 / n
    center = (phat + z ** 2 / (2 * n)) / denom
    half   = (z * math.sqrt((phat * (1 - phat) / n) + z ** 2 / (4 * n * n))) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def safe_round(x, nd=6):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    return round(float(x), nd)


# ------------------------------------------------------------------
# Performance metrics
# ------------------------------------------------------------------

def compute_metrics_with_ci(pred_R, actual_R, actual_S, universe):
    """
    Compute TP/TN/FP/FN and sensitivity/specificity/precision/NPV
    each with 95% Wilson CI.
    """
    pred_R = pred_R & universe
    pred_S = universe - pred_R

    TP = len(pred_R & actual_R)
    FP = len(pred_R & actual_S)
    TN = len(pred_S & actual_S)
    FN = len(pred_S & actual_R)

    sens = TP / (TP + FN) if (TP + FN) > 0 else None
    spec = TN / (TN + FP) if (TN + FP) > 0 else None
    prec = TP / (TP + FP) if (TP + FP) > 0 else None
    npv  = TN / (TN + FN) if (TN + FN) > 0 else None

    sens_ci = wilson_ci(TP, TP + FN) if (TP + FN) > 0 else (None, None)
    spec_ci = wilson_ci(TN, TN + FP) if (TN + FP) > 0 else (None, None)
    prec_ci = wilson_ci(TP, TP + FP) if (TP + FP) > 0 else (None, None)
    npv_ci  = wilson_ci(TN, TN + FN) if (TN + FN) > 0 else (None, None)

    return {
        "TP": TP, "TN": TN, "FP": FP, "FN": FN,
        "Sensitivity": sens, "Sens_Lower": sens_ci[0], "Sens_Upper": sens_ci[1],
        "Specificity": spec, "Spec_Lower": spec_ci[0], "Spec_Upper": spec_ci[1],
        "Precision":   prec, "Prec_Lower": prec_ci[0], "Prec_Upper": prec_ci[1],
        "NPV":         npv,  "NPV_Lower":  npv_ci[0],  "NPV_Upper":  npv_ci[1],
    }


# ------------------------------------------------------------------
# Data loading
# ------------------------------------------------------------------

def load_id_list(file):
    """
    Load ID/phenotype table. Expects columns Run and pheno.
    IDs with inconsistent phenotypes across duplicate rows are excluded.
    Only R and S phenotype values are accepted.
    """
    id_df = pd.read_csv(file, sep="\t", dtype=str)
    id_df.columns = [c.strip() for c in id_df.columns]
    id_df = id_df.rename(columns={"Run": "ID", "pheno": "Phenotype"})
    id_df.dropna(subset=["ID", "Phenotype"], inplace=True)
    id_df["ID"]        = id_df["ID"].str.strip()
    id_df["Phenotype"] = id_df["Phenotype"].str.strip().str.upper()

    invalid = set(id_df["Phenotype"].unique()) - {"R", "S"}
    if invalid:
        print(f"ERROR: Invalid phenotypes found: {invalid}. Only 'R' or 'S' are allowed.", file=sys.stderr)
        sys.exit(1)

    duplicated    = id_df[id_df.duplicated(subset=["ID"], keep=False)]
    pheno_counts  = duplicated.groupby("ID")["Phenotype"].nunique()
    inconsistent  = pheno_counts[pheno_counts > 1].index.tolist()
    if inconsistent:
        id_df = id_df[~id_df["ID"].isin(inconsistent)]
        print(f"Excluded {len(inconsistent)} IDs with inconsistent phenotypes.")

    id_df = id_df.drop_duplicates(subset=["ID"], keep="first")
    return id_df


def load_all_snp_files(id_df, snp_dir):
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
                             names=["position", "ref", "alt"], usecols=[0, 1, 2], dtype=str)
            df.dropna(subset=["position", "ref", "alt"], inplace=True)
            if df.empty:
                missing += 1
                continue
            df["position"] = df["position"].astype(str).str.strip()
            df["ref"]      = df["ref"].str.upper().str.strip()
            df["alt"]      = df["alt"].str.upper().str.strip()
            snp_dict[id_]  = set(zip(df["position"], df["ref"], df["alt"]))
        except Exception as e:
            print(f"Warning: error reading {path}: {e}")
            missing += 1

        if (i + 1) % 1000 == 0 or (i + 1) == total:
            print(f"  SNP: processed {i+1}/{total}")

    print(f"Loaded SNP (.ano) for {len(snp_dict)} IDs. Missing: {missing}.")
    return snp_dict, missing


def load_all_indel_files(id_df, indel_dir):
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
                             names=["position", "ref", "alt"], usecols=[0, 1, 2], dtype=str)
            df.dropna(subset=["position", "ref", "alt"], inplace=True)
            if df.empty:
                missing += 1
                continue
            df["position"]  = df["position"].astype(str).str.strip()
            df["ref"]       = df["ref"].str.upper().str.strip()
            df["alt"]       = df["alt"].str.upper().str.strip()
            indel_dict[id_] = set(zip(df["position"], df["ref"], df["alt"]))
        except Exception as e:
            print(f"Warning: error reading {path}: {e}")
            missing += 1

        if (i + 1) % 1000 == 0 or (i + 1) == total:
            print(f"  INDEL: processed {i+1}/{total}")

    print(f"Loaded INDEL for {len(indel_dict)} IDs. Missing: {missing}.")
    return indel_dict, missing


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)

    # Load variant list (tab-delimited with header)
    variants_df = pd.read_csv(args.variants_file, sep="\t", dtype=str)
    variants_df.columns = [c.strip() for c in variants_df.columns]
    required_vcols = {"variant", "position", "ref", "alt"}
    missing_vcols  = required_vcols - set(variants_df.columns)
    if missing_vcols:
        print(f"ERROR: Variant file missing required columns: {sorted(missing_vcols)}", file=sys.stderr)
        sys.exit(1)
    variants_df.dropna(subset=["position", "ref", "alt"], inplace=True)
    variants_df["position"] = variants_df["position"].astype(str).str.strip()
    variants_df["ref"]      = variants_df["ref"].str.upper().str.strip()
    variants_df["alt"]      = variants_df["alt"].str.upper().str.strip()

    core_vcols       = {"variant", "position", "ref", "alt"}
    extra_vcols      = [c for c in variants_df.columns if c not in core_vcols]
    print(f"Loaded {len(variants_df)} variants.")

    # Load IDs and genotype files
    id_df                  = load_id_list(args.id_list_file)
    snp_dict,   miss_snp   = load_all_snp_files(id_df, args.snp_dir)
    indel_dict, miss_indel = load_all_indel_files(id_df, args.indel_dir)

    included_ids = set(snp_dict.keys()) & set(indel_dict.keys())
    id_df        = id_df[id_df["ID"].isin(included_ids)].copy()
    print(f"IDs included (SNP ∩ INDEL): {len(included_ids)}")

    if not included_ids:
        print("ERROR: No IDs have both SNP and INDEL files. Exiting.", file=sys.stderr)
        sys.exit(1)

    universe = set(id_df["ID"])
    actual_R = set(id_df[id_df["Phenotype"] == "R"]["ID"])
    actual_S = set(id_df[id_df["Phenotype"] == "S"]["ID"])

    # Build variant -> ID mapping
    variant_to_ids = defaultdict(set)
    for id_ in universe:
        for v in snp_dict[id_] | indel_dict[id_]:
            variant_to_ids[v].add(id_)

    # Full-list predictions
    variant_keys = list(zip(variants_df["position"], variants_df["ref"], variants_df["alt"]))
    variant_set  = set(variant_keys)

    full_pred_R = set()
    for v in variant_set:
        full_pred_R |= variant_to_ids.get(v, set())

    # 1. Overall metrics
    overall    = compute_metrics_with_ci(full_pred_R, actual_R, actual_S, universe)
    overall_df = pd.DataFrame([{
        "N_isolates":  len(universe),
        "N_R":         len(actual_R),
        "N_S":         len(actual_S),
        "TP": overall["TP"], "TN": overall["TN"],
        "FP": overall["FP"], "FN": overall["FN"],
        "Sensitivity": safe_round(overall["Sensitivity"]),
        "Sens_Lower":  safe_round(overall["Sens_Lower"]),
        "Sens_Upper":  safe_round(overall["Sens_Upper"]),
        "Specificity": safe_round(overall["Specificity"]),
        "Spec_Lower":  safe_round(overall["Spec_Lower"]),
        "Spec_Upper":  safe_round(overall["Spec_Upper"]),
        "Precision":   safe_round(overall["Precision"]),
        "Prec_Lower":  safe_round(overall["Prec_Lower"]),
        "Prec_Upper":  safe_round(overall["Prec_Upper"]),
        "NPV":         safe_round(overall["NPV"]),
        "NPV_Lower":   safe_round(overall["NPV_Lower"]),
        "NPV_Upper":   safe_round(overall["NPV_Upper"]),
    }])
    overall_out = os.path.join(args.output_dir, "overall_metrics.tsv")
    overall_df.to_csv(overall_out, sep="\t", index=False)
    print(f"Overall metrics saved to: {overall_out}")

    # 2. Per-variant leave-one-out analysis
    rows       = []
    n_variants = len(variants_df)

    for i, vrow in variants_df.iterrows():
        variant_name = vrow["variant"]
        pos, ref, alt = vrow["position"], vrow["ref"], vrow["alt"]
        v = (pos, ref, alt)

        ids_with_v = variant_to_ids.get(v, set()) & universe
        r_count    = len(ids_with_v & actual_R)
        s_count    = len(ids_with_v & actual_S)

        other_coverage  = set().union(*(variant_to_ids.get(ov, set()) for ov in variant_set if ov != v))
        exclusively_by  = ids_with_v - other_coverage
        loo_pred_R      = full_pred_R - exclusively_by
        r_only_by_this  = len(exclusively_by & actual_R)

        loo = compute_metrics_with_ci(loo_pred_R, actual_R, actual_S, universe)

        def delta(new, base):
            if new is None or base is None:
                return None
            return safe_round(new - base)

        rows.append({
            "Variant":               variant_name,
            "Position":              pos,
            "Ref":                   ref,
            "Alt":                   alt,
            **{c: vrow.get(c, None) for c in extra_vcols},
            "R_count":               r_count,
            "S_count":               s_count,
            "R_only_by_this_variant": r_only_by_this,
            "LOO_Sensitivity":       safe_round(loo["Sensitivity"]),
            "LOO_Specificity":       safe_round(loo["Specificity"]),
            "LOO_Precision":         safe_round(loo["Precision"]),
            "LOO_NPV":               safe_round(loo["NPV"]),
            "Delta_Sensitivity":     delta(loo["Sensitivity"], overall["Sensitivity"]),
            "Delta_Specificity":     delta(loo["Specificity"], overall["Specificity"]),
            "Delta_Precision":       delta(loo["Precision"],   overall["Precision"]),
            "Delta_NPV":             delta(loo["NPV"],         overall["NPV"]),
        })

        if (i + 1) % 100 == 0 or (i + 1) == n_variants:
            print(f"  Leave-one-out: processed {i+1}/{n_variants} variants.")

    per_variant_df = (
        pd.DataFrame(rows)
        .sort_values("Delta_Sensitivity", ascending=True)
        .reset_index(drop=True)
    )
    per_variant_out = os.path.join(args.output_dir, "per_variant_analysis.tsv")
    per_variant_df.to_csv(per_variant_out, sep="\t", index=False)
    print(f"Per-variant analysis saved to: {per_variant_out}")

    # 3. Per-isolate predictions
    id_df_out = id_df.copy()
    id_df_out["Prediction"] = id_df_out["ID"].apply(lambda i: "R" if i in full_pred_R else "S")
    predictions_out = os.path.join(args.output_dir, "isolate_predictions.tsv")
    id_df_out.to_csv(predictions_out, sep="\t", index=False)
    print(f"Isolate predictions saved to: {predictions_out}")

    # Summary
    print(f"\n=== Overall metrics (N={len(universe)}, R={len(actual_R)}, S={len(actual_S)}) ===")
    print(f"  Sensitivity: {safe_round(overall['Sensitivity'])} "
          f"[{safe_round(overall['Sens_Lower'])}, {safe_round(overall['Sens_Upper'])}]")
    print(f"  Specificity: {safe_round(overall['Specificity'])} "
          f"[{safe_round(overall['Spec_Lower'])}, {safe_round(overall['Spec_Upper'])}]")
    print(f"  Precision:   {safe_round(overall['Precision'])} "
          f"[{safe_round(overall['Prec_Lower'])}, {safe_round(overall['Prec_Upper'])}]")
    print(f"  NPV:         {safe_round(overall['NPV'])} "
          f"[{safe_round(overall['NPV_Lower'])}, {safe_round(overall['NPV_Upper'])}]")
    print(f"\nMissing SNP files:   {miss_snp}")
    print(f"Missing INDEL files: {miss_indel}")


if __name__ == "__main__":
    main()
