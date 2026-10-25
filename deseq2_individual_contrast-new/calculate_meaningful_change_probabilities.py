#!/usr/bin/env python3
"""
Calculate Meaningful Change Probabilities on DESeq2 3LRE Meta-Analysis Results.

This script runs on top of the outputs from DESeq2_individual_analysis-new.R.
For a specified threshold delta, it computes the predictive probability that a
hypothetical new study will observe an effect size |LFC| > delta:

    P_up   = 1 - Phi((delta - mu_hat) / tau_hat)
    P_down = Phi((-delta - mu_hat) / tau_hat)
    P_any  = P_up + P_down

where:
    - mu_hat is the pooled log2 fold-change (Meta_3LRE_LFC)
    - tau_hat is the total heterogeneity standard deviation:
        tau_hat = sqrt(tau2_Level3 + tau2_Level2)
      (Level 3: between-study variance; Level 2: within-study / between-contrast variance)
    - Phi is the standard normal cumulative distribution function (CDF)
    - "up" and "down" directions follow the annotation-based sign multiplier
      (positive values indicate upregulation during biofilm formation relative to control).

The script appends these probability metrics directly into DESeq2_contrast_results.csv
for all target results directories and prints a comprehensive summary report.
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
from scipy.special import ndtr

# Ensure stdout and stderr handle UTF-8 cleanly on Windows consoles
if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass


def format_delta_col(prefix: str, delta: float) -> str:
    """Format column name with delta value (e.g. delta_0.5, delta_1.0)."""
    if delta == round(delta, 1):
        d_str = f"{delta:.1f}"
    else:
        d_str = f"{delta:g}"
    return f"{prefix}_delta_{d_str}"


def calculate_probabilities_for_df(df: pd.DataFrame, deltas: list[float]) -> tuple[pd.DataFrame, dict]:
    """
    Calculate Meta_3LRE_tau and P_up, P_down, P_any for each delta.
    Returns the updated DataFrame and a dictionary of calculation statistics.
    """
    # Verify required input columns exist
    req_cols = ["Meta_3LRE_LFC", "Heterogeneity_3L_tau2_Level3", "Heterogeneity_3L_tau2_Level2"]
    for col in req_cols:
        if col not in df.columns:
            raise KeyError(f"Required column '{col}' not found in results dataframe.")

    # Clean up any existing probability/tau columns to ensure idempotency
    existing_prob_cols = [c for c in df.columns if c == "Meta_3LRE_tau" or (c.startswith("P_") and "_delta_" in c)]
    if existing_prob_cols:
        df.drop(columns=existing_prob_cols, inplace=True)

    mu = df["Meta_3LRE_LFC"].astype(float).to_numpy()
    tau2_l3 = df["Heterogeneity_3L_tau2_Level3"].astype(float).to_numpy()
    tau2_l2 = df["Heterogeneity_3L_tau2_Level2"].astype(float).to_numpy()

    # Total heterogeneity variance: tau^2 = tau2_Level3 + tau2_Level2
    tau2_total = np.where(
        np.isnan(tau2_l3) | np.isnan(tau2_l2),
        np.nan,
        np.maximum(0.0, tau2_l3 + tau2_l2)
    )
    tau = np.sqrt(tau2_total)

    # Valid mask
    valid_mask = ~np.isnan(mu) & ~np.isnan(tau)

    # Add total tau column
    df["Meta_3LRE_tau"] = tau

    results_stats = {
        "total_genes": len(df),
        "valid_genes": int(valid_mask.sum()),
        "deltas": {}
    }

    for delta in deltas:
        delta = float(delta)
        col_up = format_delta_col("P_up", delta)
        col_down = format_delta_col("P_down", delta)
        col_any = format_delta_col("P_any", delta)

        p_up = np.full(len(df), np.nan, dtype=float)
        p_down = np.full(len(df), np.nan, dtype=float)
        p_any = np.full(len(df), np.nan, dtype=float)

        # Identify degenerate cases where tau is virtually zero (tau < 1e-12)
        degen_mask = valid_mask & (tau < 1e-12)
        norm_mask = valid_mask & (tau >= 1e-12)

        # Degenerate case handling (step function)
        if np.any(degen_mask):
            mu_deg = mu[degen_mask]
            p_up[degen_mask] = np.where(mu_deg > delta, 1.0, 0.0)
            p_down[degen_mask] = np.where(mu_deg < -delta, 1.0, 0.0)
            p_any[degen_mask] = p_up[degen_mask] + p_down[degen_mask]

        # Standard normal CDF calculation
        # P_up = 1 - Phi((delta - mu) / tau) = Phi((mu - delta) / tau)
        # P_down = Phi((-delta - mu) / tau)
        if np.any(norm_mask):
            mu_norm = mu[norm_mask]
            tau_norm = tau[norm_mask]

            z_up_arg = (mu_norm - delta) / tau_norm
            z_down_arg = (-delta - mu_norm) / tau_norm

            p_up_calc = ndtr(z_up_arg)
            p_down_calc = ndtr(z_down_arg)
            p_any_calc = p_up_calc + p_down_calc

            # Numerical safety bounds [0.0, 1.0]
            p_up[norm_mask] = np.clip(p_up_calc, 0.0, 1.0)
            p_down[norm_mask] = np.clip(p_down_calc, 0.0, 1.0)
            p_any[norm_mask] = np.clip(p_any_calc, 0.0, 1.0)

        df[col_up] = p_up
        df[col_down] = p_down
        df[col_any] = p_any

        # Compute summary stats for valid genes
        valid_p_up = p_up[valid_mask]
        valid_p_down = p_down[valid_mask]
        valid_p_any = p_any[valid_mask]

        results_stats["deltas"][delta] = {
            "col_up": col_up,
            "col_down": col_down,
            "col_any": col_any,
            "mean_p_up": float(np.mean(valid_p_up)) if len(valid_p_up) > 0 else np.nan,
            "median_p_up": float(np.median(valid_p_up)) if len(valid_p_up) > 0 else np.nan,
            "mean_p_down": float(np.mean(valid_p_down)) if len(valid_p_down) > 0 else np.nan,
            "median_p_down": float(np.median(valid_p_down)) if len(valid_p_down) > 0 else np.nan,
            "mean_p_any": float(np.mean(valid_p_any)) if len(valid_p_any) > 0 else np.nan,
            "median_p_any": float(np.median(valid_p_any)) if len(valid_p_any) > 0 else np.nan,
            "count_p_up_ge_0_8": int((valid_p_up >= 0.8).sum()),
            "count_p_up_ge_0_5": int((valid_p_up >= 0.5).sum()),
            "count_p_up_ge_0_2": int((valid_p_up >= 0.2).sum()),
            "count_p_down_ge_0_8": int((valid_p_down >= 0.8).sum()),
            "count_p_down_ge_0_5": int((valid_p_down >= 0.5).sum()),
            "count_p_down_ge_0_2": int((valid_p_down >= 0.2).sum()),
            "count_p_any_ge_0_8": int((valid_p_any >= 0.8).sum()),
            "count_p_any_ge_0_5": int((valid_p_any >= 0.5).sum()),
        }

    return df, results_stats


def generate_report_text(dir_name: str, df: pd.DataFrame, stats: dict, deltas: list[float], top_n: int = 10) -> str:
    """Generate human-readable summary report."""
    lines = []
    sep = "=" * 80
    sub_sep = "-" * 80
    lines.append(sep)
    lines.append(f"PREDICTIVE PROBABILITY OF MEANINGFUL CHANGE IN NEW STUDY (|LFC| > delta)")
    lines.append(f"Dataset Directory: {dir_name}")
    lines.append(sep)
    lines.append(f"Total evaluated clusters/genes: {stats['total_genes']}")
    lines.append(f"Genes with valid 3LRE model:   {stats['valid_genes']} ({stats['valid_genes']/stats['total_genes']*100:.1f}%)")
    lines.append(f"Sign Convention:              Positive LFC = Upregulated in Biofilm vs Control")
    lines.append(f"Heterogeneity SD formula:      tau_hat = sqrt(tau2_Level3 + tau2_Level2)")
    lines.append("")

    gene_label_col = "Representative_Gene" if "Representative_Gene" in df.columns else "cluster_id"

    for delta in deltas:
        d_stats = stats["deltas"][delta]
        col_up = d_stats["col_up"]
        col_down = d_stats["col_down"]
        col_any = d_stats["col_any"]
        n_val = stats["valid_genes"]

        lines.append(sub_sep)
        lines.append(f"Threshold delta = {delta} (|LFC| > {delta}, corresponding to {2**delta:.2f}-fold change)")
        lines.append(sub_sep)
        lines.append(f"  P_up   (Biofilm Increase): Mean = {d_stats['mean_p_up']:.4f}, Median = {d_stats['median_p_up']:.4f}")
        lines.append(f"  P_down (Biofilm Decrease): Mean = {d_stats['mean_p_down']:.4f}, Median = {d_stats['median_p_down']:.4f}")
        lines.append(f"  P_any  (Any Change):       Mean = {d_stats['mean_p_any']:.4f}, Median = {d_stats['median_p_any']:.4f}")
        lines.append("")
        lines.append("  Confidence Distribution Breakdown:")
        lines.append(f"    - High Probability Increase (P_up >= 0.80):    {d_stats['count_p_up_ge_0_8']:>5} genes ({d_stats['count_p_up_ge_0_8']/n_val*100:.2f}%)")
        lines.append(f"    - Moderate Probability Inc. (P_up >= 0.50):    {d_stats['count_p_up_ge_0_5']:>5} genes ({d_stats['count_p_up_ge_0_5']/n_val*100:.2f}%)")
        lines.append(f"    - Possible Probability Inc. (P_up >= 0.20):    {d_stats['count_p_up_ge_0_2']:>5} genes ({d_stats['count_p_up_ge_0_2']/n_val*100:.2f}%)")
        lines.append(f"    - High Probability Decrease (P_down >= 0.80):  {d_stats['count_p_down_ge_0_8']:>5} genes ({d_stats['count_p_down_ge_0_8']/n_val*100:.2f}%)")
        lines.append(f"    - Moderate Probability Dec. (P_down >= 0.50):  {d_stats['count_p_down_ge_0_5']:>5} genes ({d_stats['count_p_down_ge_0_5']/n_val*100:.2f}%)")
        lines.append(f"    - Possible Probability Dec. (P_down >= 0.20):  {d_stats['count_p_down_ge_0_2']:>5} genes ({d_stats['count_p_down_ge_0_2']/n_val*100:.2f}%)")
        lines.append(f"    - High Probability Any Change (P_any >= 0.80): {d_stats['count_p_any_ge_0_8']:>5} genes ({d_stats['count_p_any_ge_0_8']/n_val*100:.2f}%)")
        lines.append(f"    - Moderate Prob. Any Change (P_any >= 0.50):   {d_stats['count_p_any_ge_0_5']:>5} genes ({d_stats['count_p_any_ge_0_5']/n_val*100:.2f}%)")
        lines.append("")

        # Top genes for P_up
        valid_df = df[df[col_up].notna()].copy()
        if not valid_df.empty:
            top_up = valid_df.sort_values(by=col_up, ascending=False).head(top_n)
            lines.append(f"  Top {min(top_n, len(top_up))} Genes Most Likely to INCREASE in Biofilm (P_up, delta={delta}):")
            header_str = f"    {'Rank':<4} {'Gene/Cluster':<20} {'Meta_3LRE_LFC':>14} {'tau_hat':>10} {'P_up':>10} {'P_down':>10} {'P_any':>10}"
            lines.append(header_str)
            lines.append("    " + "-" * (len(header_str) - 4))
            for rank, (_, row) in enumerate(top_up.iterrows(), 1):
                gene_name = str(row.get(gene_label_col, row.get("cluster_id", "NA")))[:18]
                lines.append(
                    f"    {rank:<4} {gene_name:<20} {row['Meta_3LRE_LFC']:>14.4f} {row['Meta_3LRE_tau']:>10.4f} "
                    f"{row[col_up]:>10.4f} {row[col_down]:>10.4f} {row[col_any]:>10.4f}"
                )
            lines.append("")

            # Top genes for P_down
            top_down = valid_df.sort_values(by=col_down, ascending=False).head(top_n)
            lines.append(f"  Top {min(top_n, len(top_down))} Genes Most Likely to DECREASE in Biofilm (P_down, delta={delta}):")
            lines.append(header_str)
            lines.append("    " + "-" * (len(header_str) - 4))
            for rank, (_, row) in enumerate(top_down.iterrows(), 1):
                gene_name = str(row.get(gene_label_col, row.get("cluster_id", "NA")))[:18]
                lines.append(
                    f"    {rank:<4} {gene_name:<20} {row['Meta_3LRE_LFC']:>14.4f} {row['Meta_3LRE_tau']:>10.4f} "
                    f"{row[col_up]:>10.4f} {row[col_down]:>10.4f} {row[col_any]:>10.4f}"
                )
            lines.append("")

            # Top genes for P_any
            top_any = valid_df.sort_values(by=col_any, ascending=False).head(top_n)
            lines.append(f"  Top {min(top_n, len(top_any))} Genes Most Likely to CHANGE IN EITHER DIRECTION (|LFC| > {delta}):")
            lines.append(header_str)
            lines.append("    " + "-" * (len(header_str) - 4))
            for rank, (_, row) in enumerate(top_any.iterrows(), 1):
                gene_name = str(row.get(gene_label_col, row.get("cluster_id", "NA")))[:18]
                lines.append(
                    f"    {rank:<4} {gene_name:<20} {row['Meta_3LRE_LFC']:>14.4f} {row['Meta_3LRE_tau']:>10.4f} "
                    f"{row[col_up]:>10.4f} {row[col_down]:>10.4f} {row[col_any]:>10.4f}"
                )
            lines.append("")

    lines.append(sep)
    return "\n".join(lines)


def process_results_directory(results_dir: str, deltas: list[float], top_n: int = 10, make_backup: bool = True) -> bool:
    """Process a single results directory."""
    csv_path = os.path.join(results_dir, "DESeq2_contrast_results.csv")
    if not os.path.isfile(csv_path):
        print(f"[-] Skipping '{results_dir}': DESeq2_contrast_results.csv not found.")
        return False

    print(f"\n[*] Processing directory: {results_dir}")
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"[-] Error reading {csv_path}: {e}")
        return False

    if "Meta_3LRE_LFC" not in df.columns:
        print(f"[-] Skipping '{results_dir}': 'Meta_3LRE_LFC' column missing from {csv_path}.")
        return False

    # Create backup if requested
    if make_backup:
        bak_path = csv_path + ".bak"
        if not os.path.exists(bak_path):
            try:
                import shutil
                shutil.copy2(csv_path, bak_path)
                print(f"    Created backup: {os.path.basename(bak_path)}")
            except Exception as e:
                print(f"    Warning: Could not create backup: {e}")

    # Calculate probabilities
    updated_df, stats = calculate_probabilities_for_df(df, deltas)

    # Save updated file
    try:
        updated_df.to_csv(csv_path, index=False, encoding="utf-8")
        print(f"    Updated: {csv_path}")
    except Exception as e:
        print(f"[-] Error saving {csv_path}: {e}")
        return False

    # Also update significant 3LRE subset file if it exists and has rows
    sig_3lre_path = os.path.join(results_dir, "DESeq2_contrast_results_significant_Meta3LRE.csv")
    if os.path.isfile(sig_3lre_path):
        try:
            sig_df = pd.read_csv(sig_3lre_path)
            if not sig_df.empty and "cluster_id" in sig_df.columns:
                # Merge probability columns
                prob_cols = ["cluster_id", "Meta_3LRE_tau"] + [
                    format_delta_col(p, d) for d in deltas for p in ["P_up", "P_down", "P_any"]
                ]
                subset_probs = updated_df[[c for c in prob_cols if c in updated_df.columns]]
                # Drop existing prob cols if already in sig_df
                cols_to_drop = [c for c in prob_cols if c in sig_df.columns and c != "cluster_id"]
                if cols_to_drop:
                    sig_df.drop(columns=cols_to_drop, inplace=True)
                merged_sig = pd.merge(sig_df, subset_probs, on="cluster_id", how="left")
                merged_sig.to_csv(sig_3lre_path, index=False, encoding="utf-8")
                print(f"    Updated significant subset file: {sig_3lre_path}")
        except Exception as e:
            print(f"    Warning: Could not update significant subset file: {e}")

    # Generate and save report
    report_text = generate_report_text(os.path.basename(results_dir), updated_df, stats, deltas, top_n)
    summary_path = os.path.join(results_dir, "DESeq2_probability_metrics_summary.txt")
    try:
        with open(summary_path, "w", encoding="utf-8") as f:
            f.write(report_text)
        print(f"    Saved summary report: {summary_path}")
    except Exception as e:
        print(f"    Warning: Could not write summary file: {e}")

    # Print summary to console
    print(report_text)
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Calculate P_up, P_down, P_any predictive probabilities of meaningful change for DESeq2 meta-analysis."
    )
    parser.add_argument(
        "--results_dir",
        type=str,
        default=None,
        help="Path to specific results directory. If not specified, automatically processes all results_sample_info-* directories."
    )
    parser.add_argument(
        "--deltas",
        type=str,
        default="0.5,1.0",
        help="Comma-separated threshold delta values in log2 fold-change units (default: '0.5,1.0')."
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=10,
        help="Number of top genes to display per metric in summary reports (default: 10)."
    )
    parser.add_argument(
        "--no_backup",
        action="store_true",
        default=False,
        help="Disable automatic .bak creation before overwriting CSV files."
    )

    args = parser.parse_args()

    # Parse deltas
    try:
        deltas = [float(x.strip()) for x in args.deltas.split(",") if x.strip()]
        if not deltas:
            raise ValueError("No valid delta values found.")
    except Exception as e:
        print(f"Error parsing --deltas '{args.deltas}': {e}", file=sys.stderr)
        sys.exit(1)

    base_dir = os.path.abspath(os.path.dirname(__file__))

    # Identify directories to process
    if args.results_dir:
        target_dirs = [os.path.abspath(args.results_dir)]
    else:
        # Auto-discover results_sample_info-* or results_* containing DESeq2_contrast_results.csv
        pattern = os.path.join(base_dir, "results_*")
        candidates = sorted(glob.glob(pattern))
        target_dirs = [
            d for d in candidates
            if os.path.isdir(d) and os.path.isfile(os.path.join(d, "DESeq2_contrast_results.csv"))
        ]

    if not target_dirs:
        print("No matching results directories with DESeq2_contrast_results.csv found.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(target_dirs)} directory(ies) to process.")
    success_count = 0
    for rdir in target_dirs:
        if process_results_directory(rdir, deltas, top_n=args.top_n, make_backup=not args.no_backup):
            success_count += 1

    print(f"\n[+] Successfully completed processing {success_count}/{len(target_dirs)} directories.")


if __name__ == "__main__":
    main()
