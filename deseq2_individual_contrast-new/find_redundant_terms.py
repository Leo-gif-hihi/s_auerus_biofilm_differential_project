"""
Find redundant (overlapping) GO terms between UP and DOWN ORA significant terms.
Analyzes across all result directories containing the target cluster CSV files:
  - ORA_Meta_Analysis_Significant_Term_Clusters_up.csv
  - ORA_Meta_Analysis_Significant_Term_Clusters_down.csv
"""

import os
import sys
import argparse
from pathlib import Path
from typing import List, Optional, Dict, Any
import pandas as pd

TARGET_UP_FILE = "ORA_Meta_Analysis_Significant_Term_Clusters_up.csv"
TARGET_DOWN_FILE = "ORA_Meta_Analysis_Significant_Term_Clusters_down.csv"
TARGET_FILES = [TARGET_UP_FILE, TARGET_DOWN_FILE]


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
                f"  [!] Warning: '{file_path.name}' is locked by another program (e.g. Excel). "
                f"Saved to '{alt_path.name}' instead.",
                file=sys.stderr,
            )
            return alt_path
        except Exception as e:
            print(f"  [-] Error saving '{file_path.name}': {e}", file=sys.stderr)
            return None


def find_target_directories(search_paths: Optional[List[str]] = None) -> List[Path]:
    """
    Search for all directories containing both target ORA cluster CSV files.

    If search_paths is specified:
      - Searches each specified directory and its subdirectories.
    If search_paths is empty or None:
      - Searches recursively starting from the current working directory.
      - If none found, also searches from the parent directory.

    Skips hidden directories and python caches (__pycache__).
    """
    skip_dirs = {".git", "__pycache__", ".gemini", ".venv", "venv", ".Rproj.user"}
    matched = set()

    def scan_tree(root_path: Path):
        if not root_path.is_dir():
            return
        for root, dirs, _ in os.walk(root_path):
            dirs[:] = [d for d in dirs if d not in skip_dirs and not d.startswith(".")]
            p = Path(root)
            if all((p / tf).is_file() for tf in TARGET_FILES):
                matched.add(p.resolve())

    if search_paths:
        for sp in search_paths:
            p = Path(sp)
            if p.is_dir():
                scan_tree(p)
            else:
                print(f"[!] Warning: Specified path '{sp}' is not a directory.", file=sys.stderr)
    else:
        scan_tree(Path.cwd())
        if not matched:
            scan_tree(Path.cwd().parent)

    return sorted(list(matched))


def analyze_redundancy(data_dir: Path, export_csv: bool = True) -> Dict[str, Any]:
    """
    Analyze GO term redundancy between UP and DOWN significant terms in a single directory.
    """
    up_file = data_dir / TARGET_UP_FILE
    down_file = data_dir / TARGET_DOWN_FILE

    if not up_file.exists():
        raise FileNotFoundError(f"File not found: {up_file}")
    if not down_file.exists():
        raise FileNotFoundError(f"File not found: {down_file}")

    print("\n" + "=" * 76)
    print("ORA SIGNIFICANT TERMS REDUNDANCY ANALYSIS (UP vs DOWN)")
    print("=" * 76)
    print(f"Data directory: {data_dir.resolve()}")
    print(f"UP file:        {up_file.name}")
    print(f"DOWN file:      {down_file.name}\n")

    df_up = pd.read_csv(up_file)
    df_down = pd.read_csv(down_file)

    n_up = len(df_up)
    n_down = len(df_down)

    set_up_ids = set(df_up["ID"]) if "ID" in df_up.columns else set()
    set_down_ids = set(df_down["ID"]) if "ID" in df_down.columns else set()

    common_ids = set_up_ids.intersection(set_down_ids)
    only_up_ids = set_up_ids - set_down_ids
    only_down_ids = set_down_ids - set_up_ids
    all_unique_ids = set_up_ids.union(set_down_ids)

    n_common = len(common_ids)
    n_only_up = len(only_up_ids)
    n_only_down = len(only_down_ids)
    n_union = len(all_unique_ids)

    jaccard = (n_common / n_union) * 100 if n_union > 0 else 0.0
    overlap_coef = (
        (n_common / min(n_up, n_down)) * 100 if min(n_up, n_down) > 0 else 0.0
    )
    pct_up = (n_common / n_up * 100) if n_up > 0 else 0.0
    pct_down = (n_common / n_down * 100) if n_down > 0 else 0.0

    print("--- OVERALL SUMMARY ---")
    print(f"Total Significant Terms in UP:         {n_up}")
    print(f"Total Significant Terms in DOWN:       {n_down}")
    print(f"Redundant (Overlapping) Terms:            {n_common}")
    print(f"Terms unique to UP only:                  {n_only_up}")
    print(f"Terms unique to DOWN only:                {n_only_down}")
    print(f"Total Distinct Terms Across Both:         {n_union}")
    print(f"Overlap Percentage (UP):                  {pct_up:.2f}%")
    print(f"Overlap Percentage (DOWN):                {pct_down:.2f}%")
    print(f"Jaccard Similarity Index:                 {jaccard:.2f}%")
    print(f"Overlap Coefficient (Szymkiewicz-Simpson):{overlap_coef:.2f}%\n")

    # Breakdown by Ontology
    ont_up = set(df_up["Ontology"]) if "Ontology" in df_up.columns else set()
    ont_down = set(df_down["Ontology"]) if "Ontology" in df_down.columns else set()
    ontologies = sorted(ont_up.union(ont_down))

    if ontologies:
        print("--- BREAKDOWN BY ONTOLOGY ---")
        print(
            f"{'Ontology':<10} | {'UP':<8} | {'DOWN':<8} | {'Redundant':<10} | {'% Redundant in UP':<18} | {'% Redundant in DOWN':<18}"
        )
        print("-" * 84)

        for ont in ontologies:
            sub_up = set(df_up[df_up["Ontology"] == ont]["ID"]) if "Ontology" in df_up.columns else set()
            sub_down = set(df_down[df_down["Ontology"] == ont]["ID"]) if "Ontology" in df_down.columns else set()
            ont_common = sub_up.intersection(sub_down)
            p_up = (len(ont_common) / len(sub_up) * 100) if len(sub_up) > 0 else 0.0
            p_down = (len(ont_common) / len(sub_down) * 100) if len(sub_down) > 0 else 0.0
            print(
                f"{ont:<10} | {len(sub_up):<8} | {len(sub_down):<8} | {len(ont_common):<10} | {p_up:>17.2f}% | {p_down:>17.2f}%"
            )
        print("-" * 84 + "\n")

    # Merge redundant terms to show both UP and DOWN statistics
    if n_common > 0 and "ID" in df_up.columns and "ID" in df_down.columns:
        redundant_df = pd.merge(
            df_up,
            df_down,
            on=["ID", "Ontology", "Description"],
            suffixes=("_up", "_down"),
        )
        if "logOR_RE_padj_up" in redundant_df.columns and "logOR_RE_padj_down" in redundant_df.columns:
            redundant_df["min_padj"] = redundant_df[
                ["logOR_RE_padj_up", "logOR_RE_padj_down"]
            ].min(axis=1)
            redundant_df = redundant_df.sort_values(
                by=["Ontology", "min_padj"]
            ).drop(columns=["min_padj"])
    else:
        redundant_df = pd.DataFrame()

    print("--- TOP 10 REDUNDANT TERMS (by UP padj) ---")
    if not redundant_df.empty:
        preview_cols = ["ID", "Ontology", "Description"]
        if "logOR_RE_padj_up" in redundant_df.columns:
            preview_cols.append("logOR_RE_padj_up")
        if "logOR_RE_padj_down" in redundant_df.columns:
            preview_cols.append("logOR_RE_padj_down")
        preview_cols = [c for c in preview_cols if c in redundant_df.columns]
        print(redundant_df[preview_cols].head(10).to_string(index=False))
    else:
        print("  (No overlapping terms found between UP and DOWN)")
    print()

    if export_csv:
        output_csv = data_dir / "ORA_Meta_Analysis_Redundant_Terms_up_vs_down.csv"
        saved = safe_to_csv(redundant_df, output_csv)
        if saved:
            print(f"Exported detailed redundant terms table to:\n -> {saved.resolve()}\n")
        print("Tip: Run 'python filter_direction_specific_terms.py' to separate terms into UP-dominant, DOWN-dominant, and mixed.")

    return {
        "data_dir": data_dir,
        "n_up": n_up,
        "n_down": n_down,
        "n_common": n_common,
        "n_only_up": n_only_up,
        "n_only_down": n_only_down,
        "n_union": n_union,
        "pct_up": pct_up,
        "pct_down": pct_down,
        "jaccard": jaccard,
        "overlap_coef": overlap_coef,
        "redundant_df": redundant_df,
    }


def print_comparison_table(results: List[Dict[str, Any]]):
    """Print cross-directory comparative summary table."""
    print("\n" + "=" * 105)
    print("COMPARATIVE SUMMARY ACROSS ALL ANALYZED RESULT DIRECTORIES (UP vs DOWN REDUNDANCY)")
    print("=" * 105)
    header = (
        f"{'Directory':<45} | {'UP':>6} | {'DOWN':>6} | {'Overlap':>7} | "
        f"{'% Overlap UP':>12} | {'% Overlap DOWN':>14} | {'Jaccard':>8}"
    )
    print(header)
    print("-" * 105)
    cwd = Path.cwd().resolve()
    for res in results:
        d = res["data_dir"].resolve()
        try:
            rel_path = str(d.relative_to(cwd))
        except ValueError:
            rel_path = d.name
        print(
            f"{rel_path:<45} | {res['n_up']:>6} | {res['n_down']:>6} | {res['n_common']:>7} | "
            f"{res['pct_up']:>11.2f}% | {res['pct_down']:>13.2f}% | {res['jaccard']:>7.2f}%"
        )
    print("=" * 105 + "\n")


def analyze_all_directories(directories: List[Path], export_csv: bool = True) -> List[Dict[str, Any]]:
    """Process a list of directories and print comparative summary."""
    results = []
    total_dirs = len(directories)
    for idx, d in enumerate(directories, 1):
        print(f"\n[{idx}/{total_dirs}] Processing: {d}")
        try:
            res = analyze_redundancy(d, export_csv=export_csv)
            results.append(res)
        except Exception as e:
            print(f"[-] Error processing '{d}': {e}", file=sys.stderr)

    if len(results) > 1:
        print_comparison_table(results)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Find redundant (overlapping) GO terms between UP and DOWN ORA significant terms across all results directories containing target files."
    )
    parser.add_argument(
        "data_dirs",
        nargs="*",
        default=[],
        help="One or more target directories to analyze (default: automatically discover all directories containing target files)."
    )
    parser.add_argument(
        "--no-export",
        action="store_true",
        help="Do not export CSV output files."
    )
    args = parser.parse_args()

    target_dirs = find_target_directories(args.data_dirs if args.data_dirs else None)

    if not target_dirs:
        print(
            f"[-] No directories containing both '{TARGET_UP_FILE}' and '{TARGET_DOWN_FILE}' were found.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"[+] Found {len(target_dirs)} directory(ies) containing target files:")
    for d in target_dirs:
        print(f"    - {d}")

    analyze_all_directories(target_dirs, export_csv=not args.no_export)


if __name__ == "__main__":
    main()
