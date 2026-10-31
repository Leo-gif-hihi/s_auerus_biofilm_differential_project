#!/usr/bin/env python3
"""
Extract Fixed Effect (FE) Model Gene Universe and Significant Gene Lists.

This utility script processes DESeq2 contrast result directories (e.g. results_sample_info-*)
and extracts gene identifiers for:
  1. The "Universe" genes: all clusters having valid (non-null) Fixed Effect model
     p-values and log2 fold-change values (Meta_FE_pval & Meta_FE_LFC).
  2. The "Significant FE" genes: clusters meeting user-defined significance thresholds:
     Meta_FE_padj < padj_threshold and |Meta_FE_LFC| > lfc_threshold.

For each category, it extracts:
  - The cluster identifier (`cluster_id`)
  - All locus tags and gene symbols inside `All_Genes_Locus_Tags` (delimited by '/')

Outputs are deduplicated, alphabetically sorted, and written one identifier per line into
plain text files inside each target results folder (default: fe_universe.txt & fe_significant.txt).
"""

import os
import sys
import glob
import argparse
from typing import Set, Tuple, List, Optional
import pandas as pd

# Handle Windows console encoding gracefully
if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass


def extract_identifiers_from_df(df: pd.DataFrame) -> List[str]:
    """
    Extract cluster_id and all items from All_Genes_Locus_Tags for each row in df.
    Returns a sorted, deduplicated list of identifiers.
    """
    unique_ids: Set[str] = set()

    for _, row in df.iterrows():
        # 1. Add cluster_id if present
        cid = row.get("cluster_id")
        if pd.notnull(cid):
            cid_str = str(cid).strip()
            if cid_str:
                unique_ids.add(cid_str)

        # 2. Add all tags from All_Genes_Locus_Tags
        tags = row.get("All_Genes_Locus_Tags")
        if pd.notnull(tags):
            for tag in str(tags).split("/"):
                clean_tag = tag.strip()
                if clean_tag:
                    unique_ids.add(clean_tag)

    return sorted(unique_ids)


def process_results_directory(
    target_dir: str,
    padj_threshold: float = 0.05,
    lfc_threshold: float = 0.5,
    universe_out_name: str = "fe_universe.txt",
    sig_out_name: str = "fe_significant.txt"
) -> bool:
    """
    Process a single results directory containing DESeq2_contrast_results.csv.
    Generates universe and significant gene ID text files.
    """
    csv_path = os.path.join(target_dir, "DESeq2_contrast_results.csv")
    if not os.path.isfile(csv_path):
        print(f"[-] Skipping '{target_dir}': 'DESeq2_contrast_results.csv' not found.")
        return False

    print(f"\n========================================================")
    print(f"Processing directory: {target_dir}")
    print(f"========================================================")

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"[-] Error reading CSV file at '{csv_path}': {e}", file=sys.stderr)
        return False

    # Check required columns
    required_cols = ["cluster_id", "All_Genes_Locus_Tags", "Meta_FE_pval", "Meta_FE_padj", "Meta_FE_LFC"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        print(f"[-] Error: Missing required column(s) in {csv_path}: {missing_cols}", file=sys.stderr)
        return False

    total_rows = len(df)

    # 1. Universe filter: non-null FE p-value and log2FC
    valid_universe_mask = df["Meta_FE_pval"].notnull() & df["Meta_FE_LFC"].notnull()
    universe_df = df[valid_universe_mask]
    universe_ids = extract_identifiers_from_df(universe_df)

    # 2. Significant filter: padj < padj_threshold and |LFC| > lfc_threshold
    sig_mask = valid_universe_mask & (df["Meta_FE_padj"] < padj_threshold) & (df["Meta_FE_LFC"].abs() > lfc_threshold)
    sig_df = df[sig_mask]
    sig_ids = extract_identifiers_from_df(sig_df)

    # Output file paths
    universe_path = os.path.join(target_dir, universe_out_name)
    sig_path = os.path.join(target_dir, sig_out_name)

    # Write files
    try:
        with open(universe_path, "w", encoding="utf-8") as f:
            for gene_id in universe_ids:
                f.write(f"{gene_id}\n")

        with open(sig_path, "w", encoding="utf-8") as f:
            for gene_id in sig_ids:
                f.write(f"{gene_id}\n")
    except Exception as e:
        print(f"[-] Error writing output files: {e}", file=sys.stderr)
        return False

    # Print summary
    print(f"  Total clusters in CSV:           {total_rows}")
    print(f"  Fixed Effect valid clusters:     {len(universe_df)}")
    print(f"  [+] Universe unique IDs written: {len(universe_ids):<6} -> {universe_out_name}")
    print(f"  Fixed Effect significant clusters: {len(sig_df)} (padj < {padj_threshold}, |LFC| > {lfc_threshold})")
    print(f"  [+] Significant unique IDs:      {len(sig_ids):<6} -> {sig_out_name}")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Extract Fixed Effect (FE) universe and significant gene lists from DESeq2 contrast results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-d", "--results-dir",
        type=str,
        default=None,
        help="Target results directory. If omitted, automatically processes all results_sample_info* subdirectories."
    )
    parser.add_argument(
        "--padj",
        type=float,
        default=0.05,
        help="Adjusted p-value (Meta_FE_padj) threshold for significance."
    )
    parser.add_argument(
        "--lfc",
        type=float,
        default=0.5,
        help="Absolute log2 fold-change (|Meta_FE_LFC|) threshold for significance."
    )
    parser.add_argument(
        "--universe-out",
        type=str,
        default="fe_universe.txt",
        help="Output filename for universe gene identifiers."
    )
    parser.add_argument(
        "--sig-out",
        type=str,
        default="fe_significant.txt",
        help="Output filename for significant gene identifiers."
    )

    args = parser.parse_args()

    base_dir = os.path.abspath(os.path.dirname(__file__))

    # Identify directories to process
    if args.results_dir:
        resolved_dir = os.path.abspath(args.results_dir)
        if not os.path.isdir(resolved_dir):
            print(f"[-] Error: Specified directory '{args.results_dir}' does not exist.", file=sys.stderr)
            sys.exit(1)
        target_dirs = [resolved_dir]
    else:
        # Auto-discover results_sample_info* or results_* containing DESeq2_contrast_results.csv
        pattern = os.path.join(base_dir, "results_*")
        candidates = sorted(glob.glob(pattern))
        target_dirs = [
            d for d in candidates
            if os.path.isdir(d) and os.path.isfile(os.path.join(d, "DESeq2_contrast_results.csv"))
        ]

    if not target_dirs:
        print("[-] No matching results directories containing 'DESeq2_contrast_results.csv' found.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(target_dirs)} directory(ies) to process.")
    print(f"Criteria: FE Universe (Meta_FE_pval & Meta_FE_LFC not null)")
    print(f"          FE Significant (Meta_FE_padj < {args.padj} & |Meta_FE_LFC| > {args.lfc})")

    success_count = 0
    for target_dir in target_dirs:
        success = process_results_directory(
            target_dir=target_dir,
            padj_threshold=args.padj,
            lfc_threshold=args.lfc,
            universe_out_name=args.universe_out,
            sig_out_name=args.sig_out
        )
        if success:
            success_count += 1

    print(f"\n========================================================")
    print(f"[+] Completed processing: {success_count}/{len(target_dirs)} directories successfully.")
    print(f"========================================================\n")


if __name__ == "__main__":
    main()

