"""
Direction-Specific Filtering and Classification for ORA Meta-Analysis Terms.
Analyzes across all result directories containing the target cluster CSV files:
  - ORA_Meta_Analysis_Significant_Term_Clusters_up.csv
  - ORA_Meta_Analysis_Significant_Term_Clusters_down.csv

Classifies all significant GO terms into:
  1. UP-Specific:           Detected only in UP analysis
  2. UP-Dominant:           Detected in both, but stronger enrichment & frequency in UP
  3. DOWN-Specific:         Detected only in DOWN analysis
  4. DOWN-Dominant:         Detected in both, but stronger enrichment & frequency in DOWN
  5. Bidirectional / Mixed: Detected in both with balanced enrichment (broad activation/repression)
"""

import os
import sys
import re
import argparse
from pathlib import Path
from typing import List, Optional, Dict, Any
import numpy as np
import pandas as pd
from scipy.stats import norm

TARGET_UP_FILE = "ORA_Meta_Analysis_Significant_Term_Clusters_up.csv"
TARGET_DOWN_FILE = "ORA_Meta_Analysis_Significant_Term_Clusters_down.csv"
TARGET_FILES = [TARGET_UP_FILE, TARGET_DOWN_FILE]

CATEGORIES = [
    "UP-Specific",
    "UP-Dominant",
    "DOWN-Specific",
    "DOWN-Dominant",
    "Bidirectional / Mixed",
]


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


def parse_events(freq_str):
    """
    Extract total gene-comparison occurrences from string like:
    'hutU (11/40); fdh (10/40); hutI (10/40)' -> 11 + 10 + 10 = 31.
    """
    if pd.isna(freq_str):
        return 0
    matches = re.findall(r"\((\d+)/\d+\)", str(freq_str))
    return sum(int(m) for m in matches)


def classify_shared_term(row, z_threshold=1.0, delta_threshold=0.20, event_ratio_margin=0.05):
    """
    Classify a term that is significant in both directions.

    A term is UP-Dominant if:
      - (Z-score >= z_threshold and event_ratio_up >= 0.5 + margin) OR
      - (delta_logOR >= delta_threshold) OR
      - (Z-score > 0 and event_ratio_up >= 0.60)

    A term is DOWN-Dominant if:
      - (Z-score <= -z_threshold and event_ratio_up <= 0.5 - margin) OR
      - (delta_logOR <= -delta_threshold) OR
      - (Z-score < 0 and event_ratio_up <= 0.40)

    Otherwise, it is 'Bidirectional / Mixed'.
    """
    z = row["z_contrast"]
    d = row["delta_logOR"]
    er_up = row["event_ratio_up"]

    if (z >= z_threshold and er_up >= (0.5 + event_ratio_margin)) or \
       (d >= delta_threshold) or \
       (z > 0 and er_up >= 0.60):
        return "UP-Dominant"
    elif (z <= -z_threshold and er_up <= (0.5 - event_ratio_margin)) or \
         (d <= -delta_threshold) or \
         (z < 0 and er_up <= 0.40):
        return "DOWN-Dominant"
    else:
        return "Bidirectional / Mixed"


def run_direction_classification(data_dir: Path, z_threshold=1.0, delta_threshold=0.20) -> pd.DataFrame:
    """
    Run direction classification on a single result directory.
    """
    up_file = data_dir / TARGET_UP_FILE
    down_file = data_dir / TARGET_DOWN_FILE

    if not up_file.exists() or not down_file.exists():
        raise FileNotFoundError(f"Missing required CSVs in {data_dir}")

    df_up = pd.read_csv(up_file)
    df_down = pd.read_csv(down_file)

    set_up = set(df_up["ID"]) if "ID" in df_up.columns else set()
    set_down = set(df_down["ID"]) if "ID" in df_down.columns else set()

    shared_ids = set_up & set_down
    only_up_ids = set_up - set_down
    only_down_ids = set_down - set_up

    core_cols = [
        "ID", "Ontology", "Description", "Direction_Classification", "delta_logOR", 
        "z_contrast", "pval_contrast", "events_up", "events_down", "event_ratio_up", 
        "gene_jaccard_overlap"
    ]
    extra_cols = [
        "logOR_RE_estimate_up", "logOR_RE_SE_up", "logOR_RE_padj_up", "Gene_Count_up", "Genes_up",
        "logOR_RE_estimate_down", "logOR_RE_SE_down", "logOR_RE_padj_down", "Gene_Count_down", "Genes_down"
    ]

    # 1. Process Shared Terms
    if shared_ids and "ID" in df_up.columns and "ID" in df_down.columns:
        df_shared_up = df_up[df_up["ID"].isin(shared_ids)]
        df_shared_down = df_down[df_down["ID"].isin(shared_ids)]

        merged = pd.merge(
            df_shared_up,
            df_shared_down,
            on=["ID", "Ontology", "Description"],
            suffixes=("_up", "_down")
        )

        # Statistical contrasts
        merged["delta_logOR"] = merged["logOR_RE_estimate_up"] - merged["logOR_RE_estimate_down"]
        merged["se_contrast"] = np.sqrt(merged["logOR_RE_SE_up"]**2 + merged["logOR_RE_SE_down"]**2)
        merged["z_contrast"] = merged["delta_logOR"] / merged["se_contrast"]
        merged["pval_contrast"] = 2 * (1 - norm.cdf(np.abs(merged["z_contrast"])))
        merged["fold_ratio_OR"] = np.exp(merged["delta_logOR"])

        # Event count analysis
        merged["events_up"] = merged["Genes_Frequency_up"].apply(parse_events) if "Genes_Frequency_up" in merged.columns else 0
        merged["events_down"] = merged["Genes_Frequency_down"].apply(parse_events) if "Genes_Frequency_down" in merged.columns else 0
        total_events = merged["events_up"] + merged["events_down"]
        merged["event_ratio_up"] = np.where(total_events > 0, merged["events_up"] / total_events, 0.5)

        # Gene overlap between up and down for this term
        def calc_gene_overlap(row):
            genes_u = set(str(row["Genes_up"]).split("/"))
            genes_d = set(str(row["Genes_down"]).split("/"))
            jaccard = len(genes_u & genes_d) / len(genes_u | genes_d) if (genes_u | genes_d) else 0
            return round(jaccard, 3)

        merged["gene_jaccard_overlap"] = merged.apply(calc_gene_overlap, axis=1)

        # Classify shared terms
        merged["Direction_Classification"] = merged.apply(
            lambda r: classify_shared_term(r, z_threshold, delta_threshold),
            axis=1
        )

        shared_master = merged[core_cols + extra_cols].copy()
    else:
        merged = pd.DataFrame(columns=[
            "ID", "Ontology", "Description", "Direction_Classification", "delta_logOR",
            "z_contrast", "pval_contrast", "fold_ratio_OR", "events_up", "events_down",
            "event_ratio_up", "gene_jaccard_overlap"
        ] + extra_cols)
        shared_master = pd.DataFrame(columns=core_cols + extra_cols)

    # 2. Process UP-Only terms
    if only_up_ids and "ID" in df_up.columns:
        df_only_up = df_up[df_up["ID"].isin(only_up_ids)].copy()
        df_only_up["Direction_Classification"] = "UP-Specific"
        df_only_up["events_up"] = df_only_up["Genes_Frequency"].apply(parse_events) if "Genes_Frequency" in df_only_up.columns else 0
        df_only_up["events_down"] = 0
        df_only_up["event_ratio_up"] = 1.0
        df_only_up["delta_logOR"] = df_only_up["logOR_RE_estimate"] if "logOR_RE_estimate" in df_only_up.columns else np.nan
        df_only_up["z_contrast"] = np.nan
        df_only_up["pval_contrast"] = np.nan
        df_only_up["gene_jaccard_overlap"] = 0.0

        up_cols = [c for c in ["logOR_RE_estimate", "logOR_RE_SE", "logOR_RE_padj", "Gene_Count", "Genes"] if c in df_only_up.columns]
        up_master = df_only_up[core_cols + up_cols].rename(columns={
            "logOR_RE_estimate": "logOR_RE_estimate_up",
            "logOR_RE_SE": "logOR_RE_SE_up",
            "logOR_RE_padj": "logOR_RE_padj_up",
            "Gene_Count": "Gene_Count_up",
            "Genes": "Genes_up"
        })
        for col in ["logOR_RE_estimate_down", "logOR_RE_SE_down", "logOR_RE_padj_down", "Gene_Count_down", "Genes_down"]:
            up_master[col] = np.nan
    else:
        up_master = pd.DataFrame(columns=core_cols + extra_cols)

    # 3. Process DOWN-Only terms
    if only_down_ids and "ID" in df_down.columns:
        df_only_down = df_down[df_down["ID"].isin(only_down_ids)].copy()
        df_only_down["Direction_Classification"] = "DOWN-Specific"
        df_only_down["events_up"] = 0
        df_only_down["events_down"] = df_only_down["Genes_Frequency"].apply(parse_events) if "Genes_Frequency" in df_only_down.columns else 0
        df_only_down["event_ratio_up"] = 0.0
        df_only_down["delta_logOR"] = -df_only_down["logOR_RE_estimate"] if "logOR_RE_estimate" in df_only_down.columns else np.nan
        df_only_down["z_contrast"] = np.nan
        df_only_down["pval_contrast"] = np.nan
        df_only_down["gene_jaccard_overlap"] = 0.0

        down_cols = [c for c in ["logOR_RE_estimate", "logOR_RE_SE", "logOR_RE_padj", "Gene_Count", "Genes"] if c in df_only_down.columns]
        down_master = df_only_down[core_cols + down_cols].rename(columns={
            "logOR_RE_estimate": "logOR_RE_estimate_down",
            "logOR_RE_SE": "logOR_RE_SE_down",
            "logOR_RE_padj": "logOR_RE_padj_down",
            "Gene_Count": "Gene_Count_down",
            "Genes": "Genes_down"
        })
        for col in ["logOR_RE_estimate_up", "logOR_RE_SE_up", "logOR_RE_padj_up", "Gene_Count_up", "Genes_up"]:
            down_master[col] = np.nan
    else:
        down_master = pd.DataFrame(columns=core_cols + extra_cols)

    # 4. Build Unified Master Table of All Significant Terms
    dfs_to_concat = [df for df in [shared_master, up_master, down_master] if not df.empty]
    if dfs_to_concat:
        all_master = pd.concat(dfs_to_concat, ignore_index=True)
    else:
        all_master = pd.DataFrame(columns=core_cols + extra_cols)

    # Sort master: UP-Specific, UP-Dominant, DOWN-Specific, DOWN-Dominant, Bidirectional
    cat_order = {
        "UP-Specific": 1,
        "UP-Dominant": 2,
        "DOWN-Specific": 3,
        "DOWN-Dominant": 4,
        "Bidirectional / Mixed": 5
    }
    if not all_master.empty:
        all_master["sort_key"] = all_master["Direction_Classification"].map(cat_order).fillna(99)
        all_master = all_master.sort_values(
            by=["sort_key", "Ontology", "delta_logOR"],
            ascending=[True, True, False]
        ).drop(columns=["sort_key"])

    total_terms = len(all_master)

    # 5. Export clean datasets
    f_all = data_dir / "ORA_Significant_Terms_Direction_Classified.csv"
    f_shared = data_dir / "ORA_Significant_Terms_Shared_Classified.csv"
    f_up_clean = data_dir / "ORA_Significant_Terms_UP_Filtered.csv"
    f_down_clean = data_dir / "ORA_Significant_Terms_DOWN_Filtered.csv"
    f_mixed = data_dir / "ORA_Significant_Terms_Bidirectional_Mixed.csv"

    saved_all = safe_to_csv(all_master, f_all)
    saved_shared = safe_to_csv(merged, f_shared)

    up_filtered = all_master[all_master["Direction_Classification"].isin(["UP-Specific", "UP-Dominant"])].copy()
    down_filtered = all_master[all_master["Direction_Classification"].isin(["DOWN-Specific", "DOWN-Dominant"])].copy()
    mixed_filtered = all_master[all_master["Direction_Classification"] == "Bidirectional / Mixed"].copy()

    saved_up = safe_to_csv(up_filtered, f_up_clean)
    saved_down = safe_to_csv(down_filtered, f_down_clean)
    saved_mixed = safe_to_csv(mixed_filtered, f_mixed)

    # 6. Display Summary Console Output
    print("\n" + "=" * 80)
    print("DIRECTION-SPECIFIC CLASSIFICATION OF SIGNIFICANT GO TERMS")
    print("=" * 80)
    print(f"Data directory: {data_dir.resolve()}\n")

    vc = all_master["Direction_Classification"].value_counts() if not all_master.empty else pd.Series(dtype=int)
    summary_counts = {cat: int(vc.get(cat, 0)) for cat in CATEGORIES}

    print(f"--- OVERALL BREAKDOWN (Total {total_terms} Terms) ---")
    for cat in CATEGORIES:
        count = summary_counts[cat]
        pct = (count / total_terms * 100) if total_terms > 0 else 0.0
        print(f"  * {cat:<24} : {count:>4} terms ({pct:>5.1f}%)")

    print("\n--- COMBINED FUNCTIONAL PARTITIONS ---")
    up_total = len(up_filtered)
    down_total = len(down_filtered)
    mixed_total = len(mixed_filtered)
    pct_up_comb = (up_total / total_terms * 100) if total_terms > 0 else 0.0
    pct_down_comb = (down_total / total_terms * 100) if total_terms > 0 else 0.0
    pct_mixed_comb = (mixed_total / total_terms * 100) if total_terms > 0 else 0.0
    print(f"  -> Total UP-Driven Pathways      (UP-Specific + UP-Dominant)     : {up_total:>4} terms ({pct_up_comb:.1f}%)")
    print(f"  -> Total DOWN-Driven Pathways    (DOWN-Specific + DOWN-Dominant) : {down_total:>4} terms ({pct_down_comb:.1f}%)")
    print(f"  -> Truly Bidirectional / Mixed   (Broad Biofilm Stress/Rewiring) : {mixed_total:>4} terms ({pct_mixed_comb:.1f}%)")

    print("\n--- BREAKDOWN BY ONTOLOGY ACROSS CLASSES ---")
    if not all_master.empty and "Ontology" in all_master.columns:
        ct = pd.crosstab(all_master["Direction_Classification"], all_master["Ontology"])
        for ont in ["BP", "CC", "MF"]:
            if ont not in ct.columns:
                ct[ont] = 0
        for cat in CATEGORIES:
            if cat not in ct.index:
                ct.loc[cat] = 0
        ct = ct.loc[CATEGORIES, [c for c in ["BP", "CC", "MF"] if c in ct.columns]]
        print(ct.to_string())
    else:
        print("  (No terms available)")

    print("\n--- TOP UP-DOMINANT TERMS (Enriched in both, but strongly UP-biased) ---")
    top_up_dom = merged[merged["Direction_Classification"] == "UP-Dominant"].sort_values(by="delta_logOR", ascending=False)
    if len(top_up_dom) > 0:
        disp_cols = ["ID", "Ontology", "Description", "delta_logOR", "event_ratio_up", "logOR_RE_padj_up", "logOR_RE_padj_down"]
        disp_cols = [c for c in disp_cols if c in top_up_dom.columns]
        print(top_up_dom[disp_cols].head(6).to_string(index=False))
    else:
        print("  (None)")

    print("\n--- TOP DOWN-DOMINANT TERMS (Enriched in both, but strongly DOWN-biased) ---")
    top_down_dom = merged[merged["Direction_Classification"] == "DOWN-Dominant"].sort_values(by="delta_logOR", ascending=True)
    if len(top_down_dom) > 0:
        disp_cols = ["ID", "Ontology", "Description", "delta_logOR", "event_ratio_up", "logOR_RE_padj_up", "logOR_RE_padj_down"]
        disp_cols = [c for c in disp_cols if c in top_down_dom.columns]
        print(top_down_dom[disp_cols].head(6).to_string(index=False))
    else:
        print("  (None)")

    print("\n--- EXPORTED FILES ---")
    print(f"1. Master Annotated Table ({total_terms} terms) : {f_all.name}" + (" (Saved)" if saved_all else " (Failed/Locked)"))
    print(f"2. Clean UP Pathways ({up_total} terms)      : {f_up_clean.name}" + (" (Saved)" if saved_up else " (Failed/Locked)"))
    print(f"3. Clean DOWN Pathways ({down_total} terms)    : {f_down_clean.name}" + (" (Saved)" if saved_down else " (Failed/Locked)"))
    print(f"4. Bidirectional / Mixed ({mixed_total} terms)  : {f_mixed.name}" + (" (Saved)" if saved_mixed else " (Failed/Locked)"))
    print(f"5. Shared Terms Diagnostic ({len(merged)} terms): {f_shared.name}" + (" (Saved)" if saved_shared else " (Failed/Locked)"))
    print("=" * 80 + "\n")

    # Attach summary metrics for multi-directory reporting
    all_master.attrs["summary_metrics"] = {
        "data_dir": data_dir,
        "total": total_terms,
        "up_spec": summary_counts["UP-Specific"],
        "up_dom": summary_counts["UP-Dominant"],
        "down_spec": summary_counts["DOWN-Specific"],
        "down_dom": summary_counts["DOWN-Dominant"],
        "mixed": summary_counts["Bidirectional / Mixed"],
        "up_total": up_total,
        "down_total": down_total,
    }

    return all_master


def print_comparison_table(results: List[Dict[str, Any]]):
    """Print cross-directory comparative summary table for direction classification."""
    print("\n" + "=" * 115)
    print("MULTI-DIRECTORY DIRECTION-SPECIFIC CLASSIFICATION SUMMARY")
    print("=" * 115)
    header = (
        f"{'Directory':<40} | {'Total':>5} | {'UP-Spec':>7} | {'UP-Dom':>7} | "
        f"{'DOWN-Spec':>9} | {'DOWN-Dom':>8} | {'Mixed':>5} | {'% UP-Driven':>11} | {'% DOWN-Driven':>13}"
    )
    print(header)
    print("-" * 115)
    cwd = Path.cwd().resolve()
    for res in results:
        d = res["data_dir"].resolve()
        try:
            rel_path = str(d.relative_to(cwd))
        except ValueError:
            rel_path = d.name
        tot = res["total"]
        up_driven = res["up_spec"] + res["up_dom"]
        down_driven = res["down_spec"] + res["down_dom"]
        pct_up = (up_driven / tot * 100) if tot > 0 else 0.0
        pct_down = (down_driven / tot * 100) if tot > 0 else 0.0
        print(
            f"{rel_path:<40} | {tot:>5} | {res['up_spec']:>7} | {res['up_dom']:>7} | "
            f"{res['down_spec']:>9} | {res['down_dom']:>8} | {res['mixed']:>5} | "
            f"{pct_up:>10.1f}% | {pct_down:>12.1f}%"
        )
    print("=" * 115 + "\n")


def run_all_directories(
    directories: List[Path], z_threshold=1.0, delta_threshold=0.20
) -> Dict[Path, pd.DataFrame]:
    """Process a list of directories and print comparative summary."""
    output_dict = {}
    summary_list = []
    total_dirs = len(directories)

    for idx, d in enumerate(directories, 1):
        print(f"\n[{idx}/{total_dirs}] Processing: {d}")
        try:
            df = run_direction_classification(
                data_dir=d,
                z_threshold=z_threshold,
                delta_threshold=delta_threshold,
            )
            output_dict[d] = df
            if "summary_metrics" in df.attrs:
                summary_list.append(df.attrs["summary_metrics"])
        except Exception as e:
            print(f"[-] Error processing '{d}': {e}", file=sys.stderr)

    if len(summary_list) > 1:
        print_comparison_table(summary_list)

    return output_dict


def main():
    parser = argparse.ArgumentParser(
        description="Filter and classify ORA meta-analysis terms by direction specificity across all results directories containing target files."
    )
    parser.add_argument(
        "data_dirs",
        nargs="*",
        default=[],
        help="One or more target directories to analyze (default: automatically discover all directories containing target files)."
    )
    parser.add_argument(
        "--z-cutoff",
        type=float,
        default=1.0,
        help="Z-score cutoff for contrast (default: 1.0)"
    )
    parser.add_argument(
        "--delta-cutoff",
        type=float,
        default=0.20,
        help="Delta logOR cutoff (default: 0.20)"
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

    run_all_directories(
        target_dirs,
        z_threshold=args.z_cutoff,
        delta_threshold=args.delta_cutoff,
    )


if __name__ == "__main__":
    main()
