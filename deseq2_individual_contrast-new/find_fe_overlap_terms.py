#!/usr/bin/env python3
"""
Find Overlapping Significant Fixed Effect (FE) Genes/Terms Across Contrasts.

This script identifies and analyzes the overlap of significant Fixed Effect (FE)
genes/terms among three primary biological contrasts:
  1. Biofilm-to-Planktonic (B2P)
  2. Biofilm-to-Biofilm (B2B)
  3. Planktonic-to-Planktonic (P2P)

Significance is defined by:
  - Adjusted p-value (Meta_FE_padj) < padj_threshold (default: 0.05)
  - Log2 fold change (Meta_FE_LFC) relative to lfc_threshold (default: 0.5):
      * UP:        Meta_FE_LFC >= +lfc_threshold
      * DOWN:      Meta_FE_LFC <= -lfc_threshold
      * COMBINED:  |Meta_FE_LFC| >= lfc_threshold

Outputs three separate CSV result files extending with the threshold used:
  - fe_overlap_up_lfc<threshold>.csv
  - fe_overlap_down_lfc<threshold>.csv
  - fe_overlap_combined_lfc<threshold>.csv
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


def load_contrast_data(dir_path: Path, contrast_label: str) -> pd.DataFrame:
    """Load DESeq2_contrast_results.csv from target directory and validate columns."""
    csv_path = dir_path / "DESeq2_contrast_results.csv"
    if not csv_path.is_file():
        raise FileNotFoundError(
            f"Could not find 'DESeq2_contrast_results.csv' in {dir_path.resolve()}"
        )

    df = pd.read_csv(csv_path)
    required_cols = [
        "cluster_id",
        "Representative_Gene",
        "All_Genes_Locus_Tags",
        "Meta_FE_pval",
        "Meta_FE_padj",
        "Meta_FE_LFC",
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in {csv_path}: {missing}"
        )

    # Ensure cluster_id is a clean string
    df["cluster_id"] = df["cluster_id"].astype(str).str.strip()
    # Drop duplicates by cluster_id if any exist
    df = df.drop_duplicates(subset=["cluster_id"]).copy()
    return df


def filter_significant_clusters(
    df: pd.DataFrame,
    direction: str,
    padj_threshold: float,
    lfc_threshold: float,
) -> pd.DataFrame:
    """
    Filter clusters meeting significance criteria for the given direction.
    direction: 'up', 'down', or 'combined'
    """
    valid = (
        df["Meta_FE_pval"].notnull()
        & df["Meta_FE_padj"].notnull()
        & df["Meta_FE_LFC"].notnull()
    )

    if direction == "up":
        sig_mask = valid & (df["Meta_FE_padj"] < padj_threshold) & (df["Meta_FE_LFC"] >= lfc_threshold)
    elif direction == "down":
        sig_mask = valid & (df["Meta_FE_padj"] < padj_threshold) & (df["Meta_FE_LFC"] <= -lfc_threshold)
    elif direction == "combined":
        sig_mask = valid & (df["Meta_FE_padj"] < padj_threshold) & (df["Meta_FE_LFC"].abs() >= lfc_threshold)
    else:
        raise ValueError(f"Unknown direction '{direction}'. Choose 'up', 'down', or 'combined'.")

    return df[sig_mask].copy()


def compute_overlap_table(
    df_b2p_all: pd.DataFrame,
    df_b2b_all: pd.DataFrame,
    df_p2p_all: pd.DataFrame,
    direction: str,
    padj_threshold: float,
    lfc_threshold: float,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Build the comprehensive overlap table across B2P, B2B, and P2P for a specific direction.
    """
    # Index full tables by cluster_id for fast lookup
    idx_b2p = df_b2p_all.set_index("cluster_id")
    idx_b2b = df_b2b_all.set_index("cluster_id")
    idx_p2p = df_p2p_all.set_index("cluster_id")

    # Get significant subsets
    sig_b2p = filter_significant_clusters(df_b2p_all, direction, padj_threshold, lfc_threshold)
    sig_b2b = filter_significant_clusters(df_b2b_all, direction, padj_threshold, lfc_threshold)
    sig_p2p = filter_significant_clusters(df_p2p_all, direction, padj_threshold, lfc_threshold)

    set_b2p: Set[str] = set(sig_b2p["cluster_id"])
    set_b2b: Set[str] = set(sig_b2b["cluster_id"])
    set_p2p: Set[str] = set(sig_p2p["cluster_id"])

    all_sig_clusters = sorted(set_b2p | set_b2b | set_p2p)

    rows = []
    stats = {
        "total_union": len(all_sig_clusters),
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

    for cid in all_sig_clusters:
        in_b2p = cid in set_b2p
        in_b2b = cid in set_b2b
        in_p2p = cid in set_p2p

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

        # Resolve representative gene and locus tags across contrasts
        gene_name = ""
        locus_tags = ""
        for idx_df in (idx_b2p, idx_b2b, idx_p2p):
            if cid in idx_df.index:
                row_val = idx_df.loc[cid]
                if not gene_name and pd.notnull(row_val.get("Representative_Gene")):
                    g = str(row_val["Representative_Gene"]).strip()
                    if g:
                        gene_name = g
                if not locus_tags and pd.notnull(row_val.get("All_Genes_Locus_Tags")):
                    lt = str(row_val["All_Genes_Locus_Tags"]).strip()
                    if lt:
                        locus_tags = lt
            if gene_name and locus_tags:
                break

        # Contrast-specific values
        b2p_lfc = idx_b2p.loc[cid, "Meta_FE_LFC"] if cid in idx_b2p.index else np.nan
        b2p_padj = idx_b2p.loc[cid, "Meta_FE_padj"] if cid in idx_b2p.index else np.nan

        b2b_lfc = idx_b2b.loc[cid, "Meta_FE_LFC"] if cid in idx_b2b.index else np.nan
        b2b_padj = idx_b2b.loc[cid, "Meta_FE_padj"] if cid in idx_b2b.index else np.nan

        p2p_lfc = idx_p2p.loc[cid, "Meta_FE_LFC"] if cid in idx_p2p.index else np.nan
        p2p_padj = idx_p2p.loc[cid, "Meta_FE_padj"] if cid in idx_p2p.index else np.nan

        # Contrast-specific Heterogeneity I2 Total (Level 3 + Level 2) for significant contrasts only
        def extract_i2_total(idx_df: pd.DataFrame, is_sig: bool) -> float:
            if not is_sig or cid not in idx_df.index:
                return np.nan
            l3 = idx_df.loc[cid].get("Heterogeneity_3L_I2_Level3", np.nan)
            l2 = idx_df.loc[cid].get("Heterogeneity_3L_I2_Level2", np.nan)
            if pd.isnull(l3) and pd.isnull(l2):
                return np.nan
            val3 = 0.0 if pd.isnull(l3) else float(l3)
            val2 = 0.0 if pd.isnull(l2) else float(l2)
            return val3 + val2

        b2p_i2 = extract_i2_total(idx_b2p, in_b2p)
        b2b_i2 = extract_i2_total(idx_b2b, in_b2b)
        p2p_i2 = extract_i2_total(idx_p2p, in_p2p)

        # Sum of heterogeneity I2 Level 3 + Level 2 across overlapping comparisons (Overlap_Count >= 2)
        if overlap_count >= 2:
            sig_i2_vals = [v for v in [b2p_i2, b2b_i2, p2p_i2] if pd.notnull(v)]
            sum_i2_overlap = sum(sig_i2_vals) if sig_i2_vals else np.nan
        else:
            sum_i2_overlap = np.nan

        # Calculate minimum padj among significant contrasts for ranking
        sig_padjs = []
        if in_b2p and pd.notnull(b2p_padj):
            sig_padjs.append(b2p_padj)
        if in_b2b and pd.notnull(b2b_padj):
            sig_padjs.append(b2b_padj)
        if in_p2p and pd.notnull(p2p_padj):
            sig_padjs.append(p2p_padj)
        min_sig_padj = min(sig_padjs) if sig_padjs else 1.0

        entry = {
            "cluster_id": cid,
            "Representative_Gene": gene_name,
            "All_Genes_Locus_Tags": locus_tags,
            "Overlap_Category": overlap_cat,
            "Overlap_Count": overlap_count,
            "Significant_B2P": in_b2p,
            "Significant_B2B": in_b2b,
            "Significant_P2P": in_p2p,
            "B2P_FE_LFC": b2p_lfc,
            "B2P_FE_padj": b2p_padj,
            "B2P_I2_Total": b2p_i2,
            "B2B_FE_LFC": b2b_lfc,
            "B2B_FE_padj": b2b_padj,
            "B2B_I2_Total": b2b_i2,
            "P2P_FE_LFC": p2p_lfc,
            "P2P_FE_padj": p2p_padj,
            "P2P_I2_Total": p2p_i2,
            "Sum_I2_Overlap_Comparisons": sum_i2_overlap,
            "_min_sig_padj": min_sig_padj,
        }

        # For combined direction: classify direction pattern across significant contrasts
        if direction == "combined":
            sig_lfcs = []
            if in_b2p and pd.notnull(b2p_lfc):
                sig_lfcs.append(b2p_lfc)
            if in_b2b and pd.notnull(b2b_lfc):
                sig_lfcs.append(b2b_lfc)
            if in_p2p and pd.notnull(p2p_lfc):
                sig_lfcs.append(p2p_lfc)

            all_pos = all(x > 0 for x in sig_lfcs) if sig_lfcs else False
            all_neg = all(x < 0 for x in sig_lfcs) if sig_lfcs else False

            if all_pos:
                pattern = "Consistently_UP"
            elif all_neg:
                pattern = "Consistently_DOWN"
            else:
                pattern = "Discordant/Mixed"

            entry["Direction_Pattern"] = pattern

        rows.append(entry)

    res_df = pd.DataFrame(rows)

    if not res_df.empty:
        # Sort by Overlap_Count descending, then _min_sig_padj ascending
        res_df = res_df.sort_values(
            by=["Overlap_Count", "_min_sig_padj", "cluster_id"],
            ascending=[False, True, True],
        ).drop(columns=["_min_sig_padj"])
    else:
        # Empty dataframe fallback with schema
        cols = [
            "cluster_id",
            "Representative_Gene",
            "All_Genes_Locus_Tags",
            "Overlap_Category",
            "Overlap_Count",
            "Significant_B2P",
            "Significant_B2B",
            "Significant_P2P",
            "B2P_FE_LFC",
            "B2P_FE_padj",
            "B2P_I2_Total",
            "B2B_FE_LFC",
            "B2B_FE_padj",
            "B2B_I2_Total",
            "P2P_FE_LFC",
            "P2P_FE_padj",
            "P2P_I2_Total",
            "Sum_I2_Overlap_Comparisons",
        ]
        if direction == "combined":
            cols.append("Direction_Pattern")
        res_df = pd.DataFrame(columns=cols)

    return res_df, stats


def print_direction_summary(direction: str, stats: Dict[str, int], lfc: float, padj: float):
    """Print a clean formatted terminal summary for one direction."""
    d_label = direction.upper()
    tot = stats["total_union"]
    print(f"\n" + "=" * 78)
    print(f"  FE OVERLAP SUMMARY: {d_label} (padj < {padj}, |log2FC| >= {lfc})")
    print("=" * 78)
    print(f"  Total Distinct Significant Clusters (Union): {tot}")
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
        description="Find overlap significant Fixed Effect (FE) terms/genes across B2P, B2B, and P2P contrasts.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--lfc",
        "--log2fc",
        dest="lfc",
        type=float,
        default=0.5,
        help="Absolute log2 fold change threshold for significance (|Meta_FE_LFC| >= lfc).",
    )
    parser.add_argument(
        "--padj",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold for significance (Meta_FE_padj < padj).",
    )
    parser.add_argument(
        "--b2p-dir",
        type=str,
        default="results_sample_info-biofilm-to-planktonic",
        help="Directory containing Biofilm-to-Planktonic DESeq2 results.",
    )
    parser.add_argument(
        "--b2b-dir",
        type=str,
        default="results_sample_info-biofilm-to-biofilm",
        help="Directory containing Biofilm-to-Biofilm DESeq2 results.",
    )
    parser.add_argument(
        "--p2p-dir",
        type=str,
        default="results_sample_info-planktonic-to-planktonic",
        help="Directory containing Planktonic-to-Planktonic DESeq2 results.",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        type=str,
        default=".",
        help="Directory where overlap result CSV files will be saved.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="fe_overlap",
        help="Filename prefix for output CSV files.",
    )

    args = parser.parse_args()

    # Paths resolution
    base_dir = Path.cwd()
    b2p_path = (base_dir / args.b2p_dir).resolve()
    b2b_path = (base_dir / args.b2b_dir).resolve()
    p2p_path = (base_dir / args.p2p_dir).resolve()
    out_dir = (base_dir / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print("==============================================================================")
    print("FIXED EFFECT (FE) OVERLAP ANALYSIS ACROSS CONTRASTS")
    print("==============================================================================")
    print(f"Log2FC Threshold (|Meta_FE_LFC| >=): {args.lfc}")
    print(f"Padj Threshold   (Meta_FE_padj <):   {args.padj}")
    print(f"B2P Directory:   {b2p_path}")
    print(f"B2B Directory:   {b2b_path}")
    print(f"P2P Directory:   {p2p_path}")
    print(f"Output Directory:{out_dir}")
    print("==============================================================================")

    # Load data
    try:
        df_b2p = load_contrast_data(b2p_path, "Biofilm-to-Planktonic")
        df_b2b = load_contrast_data(b2b_path, "Biofilm-to-Biofilm")
        df_p2p = load_contrast_data(p2p_path, "Planktonic-to-Planktonic")
    except Exception as e:
        print(f"\n[-] Error loading input data: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"\n[+] Loaded input tables:")
    print(f"    - B2P: {len(df_b2p):>5} clusters")
    print(f"    - B2B: {len(df_b2b):>5} clusters")
    print(f"    - P2P: {len(df_p2p):>5} clusters")

    lfc_str = format_threshold(args.lfc)
    directions = ["up", "down", "combined"]
    saved_files = []

    for d in directions:
        overlap_df, stats = compute_overlap_table(
            df_b2p_all=df_b2p,
            df_b2b_all=df_b2b,
            df_p2p_all=df_p2p,
            direction=d,
            padj_threshold=args.padj,
            lfc_threshold=args.lfc,
        )

        print_direction_summary(d, stats, args.lfc, args.padj)

        out_name = f"{args.prefix}_{d}_lfc{lfc_str}.csv"
        out_file = out_dir / out_name
        saved = safe_to_csv(overlap_df, out_file)
        if saved:
            saved_files.append(saved)
            print(f"  [+] Saved {d.upper():<8} overlap table -> {out_file.name} ({len(overlap_df)} rows)")

        # Preview top all-3 or top overlapping genes in terminal
        if not overlap_df.empty:
            preview_top = overlap_df[overlap_df["Overlap_Count"] >= 2].head(5)
            if not preview_top.empty:
                cols_to_show = ["cluster_id", "Representative_Gene", "Overlap_Category", "B2P_FE_LFC", "B2B_FE_LFC", "P2P_FE_LFC", "Sum_I2_Overlap_Comparisons"]
                if "Direction_Pattern" in preview_top.columns:
                    cols_to_show.insert(3, "Direction_Pattern")
                print("\n    Sample Multi-Contrast Overlap Genes:")
                print("    " + preview_top[cols_to_show].to_string(index=False).replace("\n", "\n    "))

    print("\n" + "=" * 78)
    print(f"[+] All 3 overlap files successfully generated with threshold lfc{lfc_str}:")
    for f in saved_files:
        print(f"    * {f.name}")
    print("=" * 78 + "\n")


if __name__ == "__main__":
    main()

