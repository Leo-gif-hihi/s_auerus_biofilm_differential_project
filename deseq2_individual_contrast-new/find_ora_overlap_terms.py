#!/usr/bin/env python3
"""
Find Overlapping Significant Over-Representation Analysis (ORA) Terms Across Contrasts.

This script identifies and analyzes the overlap of significant Over-Representation
Analysis (ORA / Gene Ontology) terms among three primary biological contrasts:
  1. Biofilm-to-Planktonic (B2P)
  2. Biofilm-to-Biofilm (B2B)
  3. Planktonic-to-Planktonic (P2P)

Input files evaluated:
  - ORA_Meta_Analysis_Results_up.csv
  - ORA_Meta_Analysis_Results_down.csv
  - ORA_Meta_Analysis_Results_all.csv

Significance filtering criteria (default):
  - Model p-adjusted value (padj) < padj_threshold (default: 0.05)
  - Log Odds Ratio estimate (logOR_estimate) >= logor_threshold (default: 0.5)
  - Model choice via --model: 're' (Random-Effects, default) or 'fe' (Fixed-Effect)

Outputs three separate CSV result files extending with the threshold used:
  - ora_overlap_up_logor<threshold>.csv
  - ora_overlap_down_logor<threshold>.csv
  - ora_overlap_all_logor<threshold>.csv
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import numpy as np
import pandas as pd

# Handle Windows console encoding gracefully
if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


def safe_to_csv(df: pd.DataFrame, file_path: Path, **kwargs) -> Optional[Path]:
    """Save dataframe to CSV, handling file lock / PermissionError gracefully."""
    try:
        df.to_csv(file_path, index=False, **kwargs)
        return file_path
    except PermissionError:
        alt_path = file_path.with_name(f"{file_path.stem}_new{file_path.suffix}")
        try:
            df.to_csv(alt_path, index=False, **kwargs)
            print(
                f"  [!] Warning: '{file_path.name}' is locked by another program. "
                f"Saved to '{alt_path.name}' instead.",
                file=sys.stderr,
            )
            return alt_path
        except Exception as e:
            print(f"  [-] Error saving '{file_path.name}': {e}", file=sys.stderr)
            return None


def format_threshold(val: float) -> str:
    """Format threshold cleanly for file naming (e.g., 0.5 -> '0.5', 1.0 -> '1.0')."""
    s = f"{val:.4f}".rstrip("0")
    if s.endswith("."):
        s += "0"
    return s


def load_ora_data(dir_path: Path, direction: str, model: str) -> Optional[pd.DataFrame]:
    """
    Load ORA_Meta_Analysis_Results_{direction}.csv from target directory.
    Returns None if file is missing (e.g. down-regulated ORA in biofilm-to-planktonic).
    """
    csv_path = dir_path / f"ORA_Meta_Analysis_Results_{direction}.csv"
    if not csv_path.is_file():
        return None

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"  [-] Error reading '{csv_path}': {e}", file=sys.stderr)
        return None

    est_col = f"logOR_{model.upper()}_estimate"
    padj_col = f"logOR_{model.upper()}_padj"

    required_cols = ["ID", "Ontology", "Description", est_col, padj_col]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in {csv_path}: {missing}"
        )

    # Ensure ID is clean string
    df["ID"] = df["ID"].astype(str).str.strip()
    df = df.drop_duplicates(subset=["ID"]).copy()
    return df


def filter_significant_ora(
    df: Optional[pd.DataFrame],
    model: str,
    padj_threshold: float,
    logor_threshold: float,
) -> pd.DataFrame:
    """Filter ORA terms by model-specific padj and logOR estimate threshold."""
    if df is None or df.empty:
        return pd.DataFrame()

    est_col = f"logOR_{model.upper()}_estimate"
    padj_col = f"logOR_{model.upper()}_padj"

    valid = df[est_col].notnull() & df[padj_col].notnull()
    sig_mask = valid & (df[padj_col] < padj_threshold) & (df[est_col] >= logor_threshold)
    return df[sig_mask].copy()


def compute_ora_overlap_table(
    df_b2p_all: Optional[pd.DataFrame],
    df_b2b_all: Optional[pd.DataFrame],
    df_p2p_all: Optional[pd.DataFrame],
    model: str,
    padj_threshold: float,
    logor_threshold: float,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Build the comprehensive ORA overlap table across B2P, B2B, and P2P for a specific direction.
    """
    est_col = f"logOR_{model.upper()}_estimate"
    padj_col = f"logOR_{model.upper()}_padj"

    # Index full tables by ID for fast metric retrieval
    idx_b2p = df_b2p_all.set_index("ID") if df_b2p_all is not None else None
    idx_b2b = df_b2b_all.set_index("ID") if df_b2b_all is not None else None
    idx_p2p = df_p2p_all.set_index("ID") if df_p2p_all is not None else None

    # Get significant subsets
    sig_b2p = filter_significant_ora(df_b2p_all, model, padj_threshold, logor_threshold)
    sig_b2b = filter_significant_ora(df_b2b_all, model, padj_threshold, logor_threshold)
    sig_p2p = filter_significant_ora(df_p2p_all, model, padj_threshold, logor_threshold)

    set_b2p: Set[str] = set(sig_b2p["ID"]) if not sig_b2p.empty else set()
    set_b2b: Set[str] = set(sig_b2b["ID"]) if not sig_b2b.empty else set()
    set_p2p: Set[str] = set(sig_p2p["ID"]) if not sig_p2p.empty else set()

    all_sig_terms = sorted(set_b2p | set_b2b | set_p2p)

    stats = {
        "total_union": len(all_sig_terms),
        "b2p_total": len(set_b2p),
        "b2b_total": len(set_b2b),
        "p2p_total": len(set_p2p),
        "all_3": len(set_b2p & set_b2b & set_p2p),
        "b2p_and_b2b": len((set_b2p & set_b2b) - set_p2p),
        "b2p_and_p2p": len((set_b2p & set_p2p) - set_b2b),
        "b2b_and_p2p": len((set_b2b & set_p2p) - set_b2p),
        "b2p_only": len(set_b2p - set_b2b - set_p2p),
        "b2b_only": len(set_b2b - set_b2p - set_p2p),
        "p2p_only": len(set_p2p - set_b2p - set_b2b),
    }

    rows = []
    for term_id in all_sig_terms:
        in_b2p = term_id in set_b2p
        in_b2b = term_id in set_b2b
        in_p2p = term_id in set_p2p

        overlap_count = int(in_b2p) + int(in_b2b) + int(in_p2p)

        if overlap_count == 3:
            overlap_cat = "All_3"
        elif overlap_count == 2:
            if in_b2p and in_b2b:
                overlap_cat = "B2P_and_B2B"
            elif in_b2p and in_p2p:
                overlap_cat = "B2P_and_P2P"
            else:
                overlap_cat = "B2B_and_P2P"
        else:
            if in_b2p:
                overlap_cat = "B2P_only"
            elif in_b2b:
                overlap_cat = "B2B_only"
            else:
                overlap_cat = "P2P_only"

        # Resolve ontology and description across datasets
        ontology = ""
        description = ""
        for idx_df in (idx_b2p, idx_b2b, idx_p2p):
            if idx_df is not None and term_id in idx_df.index:
                row_val = idx_df.loc[term_id]
                if not ontology and pd.notnull(row_val.get("Ontology")):
                    ontology = str(row_val["Ontology"]).strip()
                if not description and pd.notnull(row_val.get("Description")):
                    description = str(row_val["Description"]).strip()
            if ontology and description:
                break

        # Contrast-specific values
        def get_contrast_vals(idx_df: Optional[pd.DataFrame]):
            if idx_df is not None and term_id in idx_df.index:
                r = idx_df.loc[term_id]
                est = r.get(est_col, np.nan)
                padj = r.get(padj_col, np.nan)
                count = r.get("Gene_Count", np.nan)
                genes = r.get("Genes", "")
                return est, padj, count, (genes if pd.notnull(genes) else "")
            return np.nan, np.nan, np.nan, ""

        b2p_est, b2p_padj, b2p_cnt, b2p_genes = get_contrast_vals(idx_b2p)
        b2b_est, b2b_padj, b2b_cnt, b2b_genes = get_contrast_vals(idx_b2b)
        p2p_est, p2p_padj, p2p_cnt, p2p_genes = get_contrast_vals(idx_p2p)

        # Min padj among significant contrasts for ranking
        sig_padjs = []
        if in_b2p and pd.notnull(b2p_padj):
            sig_padjs.append(b2p_padj)
        if in_b2b and pd.notnull(b2b_padj):
            sig_padjs.append(b2b_padj)
        if in_p2p and pd.notnull(p2p_padj):
            sig_padjs.append(p2p_padj)
        min_sig_padj = min(sig_padjs) if sig_padjs else 1.0

        entry = {
            "ID": term_id,
            "Ontology": ontology,
            "Description": description,
            "Overlap_Category": overlap_cat,
            "Overlap_Count": overlap_count,
            "Significant_B2P": in_b2p,
            "Significant_B2B": in_b2b,
            "Significant_P2P": in_p2p,
            f"B2P_logOR_{model.upper()}_estimate": b2p_est,
            f"B2P_logOR_{model.upper()}_padj": b2p_padj,
            "B2P_Gene_Count": b2p_cnt,
            "B2P_Genes": b2p_genes,
            f"B2B_logOR_{model.upper()}_estimate": b2b_est,
            f"B2B_logOR_{model.upper()}_padj": b2b_padj,
            "B2B_Gene_Count": b2b_cnt,
            "B2B_Genes": b2b_genes,
            f"P2P_logOR_{model.upper()}_estimate": p2p_est,
            f"P2P_logOR_{model.upper()}_padj": p2p_padj,
            "P2P_Gene_Count": p2p_cnt,
            "P2P_Genes": p2p_genes,
            "_min_sig_padj": min_sig_padj,
        }
        rows.append(entry)

    res_df = pd.DataFrame(rows)

    if not res_df.empty:
        # Sort by Overlap_Count descending, then _min_sig_padj ascending, then ID
        res_df = res_df.sort_values(
            by=["Overlap_Count", "_min_sig_padj", "ID"],
            ascending=[False, True, True],
        ).drop(columns=["_min_sig_padj"])
    else:
        cols = [
            "ID",
            "Ontology",
            "Description",
            "Overlap_Category",
            "Overlap_Count",
            "Significant_B2P",
            "Significant_B2B",
            "Significant_P2P",
            f"B2P_logOR_{model.upper()}_estimate",
            f"B2P_logOR_{model.upper()}_padj",
            "B2P_Gene_Count",
            "B2P_Genes",
            f"B2B_logOR_{model.upper()}_estimate",
            f"B2B_logOR_{model.upper()}_padj",
            "B2B_Gene_Count",
            "B2B_Genes",
            f"P2P_logOR_{model.upper()}_estimate",
            f"P2P_logOR_{model.upper()}_padj",
            "P2P_Gene_Count",
            "P2P_Genes",
        ]
        res_df = pd.DataFrame(columns=cols)

    return res_df, stats


def print_direction_summary(direction: str, stats: Dict[str, int], logor: float, padj: float, model: str):
    """Print clean formatted terminal summary for ORA overlap."""
    d_label = direction.upper()
    tot = stats["total_union"]
    print(f"\n" + "=" * 78)
    print(f"  ORA OVERLAP SUMMARY: {d_label} [{model.upper()} Model] (padj < {padj}, logOR >= {logor})")
    print("=" * 78)
    print(f"  Total Distinct Significant Terms (Union):     {tot}")
    print(f"  - Biofilm-to-Planktonic (B2P) total:         {stats['b2p_total']:>5}")
    print(f"  - Biofilm-to-Biofilm    (B2B) total:         {stats['b2b_total']:>5}")
    print(f"  - Planktonic-to-Planktonic (P2P) total:      {stats['p2p_total']:>5}")
    print("-" * 78)
    print("  Overlap Categorization Breakdown:")
    pct_all3 = (stats['all_3'] / tot * 100) if tot > 0 else 0.0
    print(f"    * Overlap across ALL 3 contrasts (All_3):  {stats['all_3']:>5} ({pct_all3:>5.1f}%)")
    print(f"    * Shared B2P & B2B only:                  {stats['b2p_and_b2b']:>5}")
    print(f"    * Shared B2P & P2P only:                  {stats['b2p_and_p2p']:>5}")
    print(f"    * Shared B2B & P2P only:                  {stats['b2b_and_p2p']:>5}")
    print(f"    * Unique to B2P only:                     {stats['b2p_only']:>5}")
    print(f"    * Unique to B2B only:                     {stats['b2b_only']:>5}")
    print(f"    * Unique to P2P only:                     {stats['p2p_only']:>5}")
    print("-" * 78)


def main():
    parser = argparse.ArgumentParser(
        description="Find overlap significant ORA terms across B2P, B2B, and P2P contrasts using logOR estimate and padj.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--logor",
        "--lfc",
        dest="logor",
        type=float,
        default=0.5,
        help="Minimum log Odds Ratio estimate threshold for significance (logOR_estimate >= logor).",
    )
    parser.add_argument(
        "--padj",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold for significance (padj < padj_threshold).",
    )
    parser.add_argument(
        "--model",
        type=str,
        default="re",
        choices=["re", "fe"],
        help="Meta-analysis model to evaluate: 're' (Random-Effects) or 'fe' (Fixed-Effect).",
    )
    parser.add_argument(
        "--b2p-dir",
        type=str,
        default="results_sample_info-biofilm-to-planktonic",
        help="Directory containing Biofilm-to-Planktonic ORA results.",
    )
    parser.add_argument(
        "--b2b-dir",
        type=str,
        default="results_sample_info-biofilm-to-biofilm",
        help="Directory containing Biofilm-to-Biofilm ORA results.",
    )
    parser.add_argument(
        "--p2p-dir",
        type=str,
        default="results_sample_info-planktonic-to-planktonic",
        help="Directory containing Planktonic-to-Planktonic ORA results.",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        type=str,
        default=".",
        help="Directory where ORA overlap result CSV files will be saved.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="ora_overlap",
        help="Filename prefix for output CSV files.",
    )

    args = parser.parse_args()

    base_dir = Path.cwd()
    b2p_path = (base_dir / args.b2p_dir).resolve()
    b2b_path = (base_dir / args.b2b_dir).resolve()
    p2p_path = (base_dir / args.p2p_dir).resolve()
    out_dir = (base_dir / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print("==============================================================================")
    print("OVER-REPRESENTATION ANALYSIS (ORA) OVERLAP ACROSS CONTRASTS")
    print("==============================================================================")
    print(f"LogOR Threshold (logOR_{args.model.upper()}_estimate >=): {args.logor}")
    print(f"Padj Threshold  (logOR_{args.model.upper()}_padj <):       {args.padj}")
    print(f"Model Selection:                                  {args.model.upper()}")
    print(f"B2P Directory:   {b2p_path}")
    print(f"B2B Directory:   {b2b_path}")
    print(f"P2P Directory:   {p2p_path}")
    print(f"Output Directory:{out_dir}")
    print("==============================================================================")

    logor_str = format_threshold(args.logor)
    directions = ["up", "down", "all"]
    saved_files = []

    for d in directions:
        # Load tables for direction
        df_b2p = load_ora_data(b2p_path, d, args.model)
        df_b2b = load_ora_data(b2b_path, d, args.model)
        df_p2p = load_ora_data(p2p_path, d, args.model)

        # Informative notice if any contrast is missing (e.g. down in B2P)
        for label, df_ref, dir_ref in [("B2P", df_b2p, b2p_path), ("B2B", df_b2b, b2b_path), ("P2P", df_p2p, p2p_path)]:
            if df_ref is None:
                print(f"  [i] Note: 'ORA_Meta_Analysis_Results_{d}.csv' not found in {label} ({dir_ref.name}); treated as 0 terms.")

        overlap_df, stats = compute_ora_overlap_table(
            df_b2p_all=df_b2p,
            df_b2b_all=df_b2b,
            df_p2p_all=df_p2p,
            model=args.model,
            padj_threshold=args.padj,
            logor_threshold=args.logor,
        )

        print_direction_summary(d, stats, args.logor, args.padj, args.model)

        out_name = f"{args.prefix}_{d}_logor{logor_str}.csv"
        out_file = out_dir / out_name
        saved = safe_to_csv(overlap_df, out_file)
        if saved:
            saved_files.append(saved)
            print(f"  [+] Saved {d.upper():<8} overlap table -> {out_file.name} ({len(overlap_df)} rows)")

        # Preview top all-3 or top overlapping terms in terminal
        if not overlap_df.empty:
            preview_top = overlap_df[overlap_df["Overlap_Count"] >= 2].head(5)
            if not preview_top.empty:
                est_c = f"logOR_{args.model.upper()}_estimate"
                cols_to_show = ["ID", "Ontology", "Description", "Overlap_Category"]
                print("\n    Sample Multi-Contrast Overlap Terms:")
                print("    " + preview_top[cols_to_show].to_string(index=False).replace("\n", "\n    "))

    print("\n" + "=" * 78)
    print(f"[+] All 3 ORA overlap files successfully generated with threshold logor{logor_str} ({args.model.upper()} model):")
    for f in saved_files:
        print(f"    * {f.name}")
    print("=" * 78 + "\n")


if __name__ == "__main__":
    main()
