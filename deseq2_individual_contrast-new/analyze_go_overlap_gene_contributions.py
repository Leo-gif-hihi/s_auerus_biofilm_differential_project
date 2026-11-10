#!/usr/bin/env python3
"""
Analyze Contribution of Overlapping Genes to Overlapping GO Terms Across Contrasts.

This script determines whether overlapping differentially expressed Fixed Effect (FE)
genes (defined as present/significant in at least 2 contrasts, e.g. Overlap_Count >= 2)
contribute to the overlapping Over-Representation Analysis (ORA) Gene Ontology (GO) terms.

Input files:
  1. GAF Annotation File:
     --gaf-file (default: QuickGO-annotations-1764241851574-20251127.gaf)
  2. FE Overlap Genes CSVs:
     --fe-pattern (default: fe_overlap_{direction}_lfc{lfc}.csv)
  3. ORA Overlap Terms CSVs:
     --ora-pattern (default: ora_overlap_{direction}_logor{logor}.csv)

Analyzes separately across three directions:
  - up:       fe_overlap_up_lfc<lfc>.csv vs ora_overlap_up_logor<logor>.csv
  - down:     fe_overlap_down_lfc<lfc>.csv vs ora_overlap_down_logor<logor>.csv
  - combined: fe_overlap_combined_lfc<lfc>.csv vs ora_overlap_all_logor<logor>.csv

Outputs three CSV result files extending with the thresholds used:
  - go_term_overlap_gene_contribution_up_logor<logor>_lfc<lfc>.csv
  - go_term_overlap_gene_contribution_down_logor<logor>_lfc<lfc>.csv
  - go_term_overlap_gene_contribution_combined_logor<logor>_lfc<lfc>.csv
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


def parse_gaf_annotations(gaf_path: Path) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
    """
    Parse GAF annotation file into a mapping:
      term2genes: GO_ID -> set of gene symbols / locus tags
      gene2name:  gene_symbol -> human-readable name (if present)
    Filters out 'NOT' qualifiers and uninformative 'ND' evidence codes.
    """
    if not gaf_path.is_file():
        raise FileNotFoundError(f"GAF file not found: {gaf_path.resolve()}")

    term2genes: Dict[str, Set[str]] = {}
    gene2name: Dict[str, str] = {}

    with open(gaf_path, "r", encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            if line.startswith("!"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 5:
                continue

            # Column 3: Qualifier (e.g. enables, NOT|enables)
            qualifier = parts[3].strip()
            if "NOT" in qualifier:
                continue

            # Column 6: Evidence code (V7 in 1-based indexing)
            evidence = parts[6].strip() if len(parts) > 6 else ""
            if evidence == "ND":
                continue

            term_id = parts[4].strip()
            symbol = parts[2].strip()
            name = parts[9].strip() if len(parts) > 9 else ""
            synonyms = (
                [s.strip() for s in parts[10].split("|") if s.strip()]
                if len(parts) > 10
                else []
            )

            if term_id not in term2genes:
                term2genes[term_id] = set()

            if symbol:
                term2genes[term_id].add(symbol)
                if name and symbol not in gene2name:
                    gene2name[symbol] = name

            for syn in synonyms:
                term2genes[term_id].add(syn)
                if name and syn not in gene2name:
                    gene2name[syn] = name

    return term2genes, gene2name


def load_overlap_genes(
    fe_file: Path, min_overlap_count: int = 2
) -> Tuple[pd.DataFrame, Dict[str, str], Dict[str, str]]:
    """
    Load FE overlap genes CSV and filter by Overlap_Count >= min_overlap_count.
    Returns:
      sig_df: filtered DataFrame of overlap genes
      tag_to_cluster: mapping of each locus tag/symbol to cluster_id
      cluster_to_label: mapping of cluster_id to display label (Gene_symbol or cluster_id)
    """
    if not fe_file.is_file():
        raise FileNotFoundError(f"FE overlap file not found: {fe_file.resolve()}")

    df = pd.read_csv(fe_file)
    required = ["cluster_id", "Overlap_Count"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {fe_file}: {missing}")

    sig_df = df[df["Overlap_Count"] >= min_overlap_count].copy()

    tag_to_cluster: Dict[str, str] = {}
    cluster_to_label: Dict[str, str] = {}

    for _, row in sig_df.iterrows():
        cid = str(row["cluster_id"]).strip()
        rep_gene = str(row["Representative_Gene"]).strip() if pd.notnull(row.get("Representative_Gene")) else ""
        label = rep_gene if rep_gene else cid
        cluster_to_label[cid] = label

        # Map cluster_id itself
        tag_to_cluster[cid] = cid

        # Map representative gene
        if rep_gene:
            tag_to_cluster[rep_gene] = cid

        # Map all locus tags
        tags_str = row.get("All_Genes_Locus_Tags")
        if pd.notnull(tags_str):
            for t in str(tags_str).split("/"):
                clean_t = t.strip()
                if clean_t:
                    tag_to_cluster[clean_t] = cid

    return sig_df, tag_to_cluster, cluster_to_label


def load_overlap_go_terms(
    ora_file: Path, min_term_overlap: int = 2, all_terms: bool = False
) -> pd.DataFrame:
    """Load ORA overlap terms CSV and optionally filter by Overlap_Count >= min_term_overlap."""
    if not ora_file.is_file():
        raise FileNotFoundError(f"ORA overlap file not found: {ora_file.resolve()}")

    df = pd.read_csv(ora_file)
    required = ["ID", "Ontology", "Description", "Overlap_Category", "Overlap_Count"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {ora_file}: {missing}")

    if not all_terms:
        df = df[df["Overlap_Count"] >= min_term_overlap].copy()

    return df


def analyze_direction(
    direction: str,
    fe_file: Path,
    ora_file: Path,
    term2genes_gaf: Dict[str, Set[str]],
    min_gene_overlap: int = 2,
    min_term_overlap: int = 2,
    all_terms: bool = False,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Analyze overlap gene contribution for a single direction."""
    try:
        genes_df, tag_to_cluster, cluster_to_label = load_overlap_genes(fe_file, min_gene_overlap)
    except FileNotFoundError as e:
        print(f"  [-] Warning: {e}; returning empty result for {direction}.", file=sys.stderr)
        return pd.DataFrame(), {}

    try:
        ora_df = load_overlap_go_terms(ora_file, min_term_overlap, all_terms)
    except FileNotFoundError as e:
        print(f"  [-] Warning: {e}; returning empty result for {direction}.", file=sys.stderr)
        return pd.DataFrame(), {}

    rows = []
    terms_with_overlap = 0

    for _, row in ora_df.iterrows():
        term_id = str(row["ID"]).strip()
        ontology = row.get("Ontology", "")
        desc = row.get("Description", "")
        term_cat = row.get("Overlap_Category", "")
        term_count = row.get("Overlap_Count", 0)

        # 1. GAF mapped genes for this term
        annotated_genes = term2genes_gaf.get(term_id, set())
        tot_annotated = len(annotated_genes)

        # 2. Identify overlap genes present in this GO term via GAF
        matched_cids: Set[str] = set()
        matched_tags: Set[str] = set()
        for gene_tag in annotated_genes:
            if gene_tag in tag_to_cluster:
                cid = tag_to_cluster[gene_tag]
                matched_cids.add(cid)
                matched_tags.add(gene_tag)

        # Sort matched overlap genes by display label
        sorted_cids = sorted(matched_cids, key=lambda c: cluster_to_label.get(c, c).lower())
        overlap_genes_str = "/".join(
            f"{cluster_to_label.get(c, c)} ({c})" if cluster_to_label.get(c, c) != c else c
            for c in sorted_cids
        )

        # 3. Active contributing DE genes from ORA results (if available in columns)
        active_de_cids: Set[str] = set()
        for col in ["B2P_Genes", "B2B_Genes", "P2P_Genes"]:
            val = row.get(col)
            if pd.notnull(val):
                for g in str(val).split("/"):
                    clean_g = g.strip()
                    if clean_g in tag_to_cluster:
                        active_de_cids.add(tag_to_cluster[clean_g])

        sorted_active_cids = sorted(active_de_cids, key=lambda c: cluster_to_label.get(c, c).lower())
        active_de_str = "/".join(
            f"{cluster_to_label.get(c, c)} ({c})" if cluster_to_label.get(c, c) != c else c
            for c in sorted_active_cids
        )

        n_contributing = len(matched_cids)
        n_active_de = len(active_de_cids)

        # Has_Overlap_Gene is True if overlap genes exist in either direct GAF annotations OR active ORA DE contributors
        has_overlap = (n_contributing > 0) or (n_active_de > 0)
        if has_overlap:
            terms_with_overlap += 1

        pct_overlap = (
            round((n_contributing / tot_annotated) * 100, 2)
            if tot_annotated > 0
            else 0.0
        )

        # Combined all-sources overlap gene list (GAF direct + ORA active DE)
        all_contributing_cids = sorted(matched_cids | active_de_cids, key=lambda c: cluster_to_label.get(c, c).lower())
        all_overlap_genes_str = "/".join(
            f"{cluster_to_label.get(c, c)} ({c})" if cluster_to_label.get(c, c) != c else c
            for c in all_contributing_cids
        )

        # 4. Compile entry
        entry = {
            "ID": term_id,
            "Ontology": ontology,
            "Description": desc,
            "Term_Overlap_Category": term_cat,
            "Term_Overlap_Count": term_count,
            "Total_Annotated_Genes": tot_annotated,
            "Has_Overlap_Gene": has_overlap,
            "All_Overlap_Genes": all_overlap_genes_str,
            "N_All_Overlap_Genes": len(all_contributing_cids),
            "N_Overlap_Genes_GAF_Direct": n_contributing,
            "Overlap_Genes_GAF_Direct": overlap_genes_str,
            "Overlap_Gene_Percentage_GAF": pct_overlap,
            "N_Active_DE_Overlap_Genes": n_active_de,
            "Active_DE_Overlap_Genes": active_de_str,
        }

        # Include contrast logOR estimates & padj if present in ORA table
        for c_prefix in ["B2P", "B2B", "P2P"]:
            for m_suffix in ["estimate", "padj"]:
                col_match = [c for c in row.index if c.startswith(f"{c_prefix}_logOR_") and c.endswith(m_suffix)]
                if col_match:
                    entry[f"{c_prefix}_logOR_{m_suffix}"] = row.get(col_match[0], np.nan)

        rows.append(entry)

    res_df = pd.DataFrame(rows)

    if not res_df.empty:
        # Sort by Has_Overlap_Gene descending, N_All_Overlap_Genes descending, Overlap_Gene_Percentage_GAF descending
        res_df = res_df.sort_values(
            by=["Has_Overlap_Gene", "N_All_Overlap_Genes", "Overlap_Gene_Percentage_GAF", "Term_Overlap_Count"],
            ascending=[False, False, False, False],
        )

    stats = {
        "total_terms_evaluated": len(ora_df),
        "total_overlap_genes_available": len(genes_df),
        "terms_with_overlap_genes": terms_with_overlap,
        "terms_without_overlap_genes": len(ora_df) - terms_with_overlap,
    }

    return res_df, stats


def print_contribution_summary(
    direction: str, stats: Dict[str, int], lfc: float, logor: float
):
    """Print terminal summary of overlap gene contribution."""
    d_label = direction.upper()
    tot_terms = stats.get("total_terms_evaluated", 0)
    with_genes = stats.get("terms_with_overlap_genes", 0)
    pct = (with_genes / tot_terms * 100) if tot_terms > 0 else 0.0

    print("\n" + "=" * 78)
    print(f"  GO OVERLAP GENE CONTRIBUTION: {d_label} (logOR >= {logor}, gene |log2FC| >= {lfc})")
    print("=" * 78)
    print(f"  Total Overlap Genes Available (Overlap_Count >= 2): {stats.get('total_overlap_genes_available', 0):>5}")
    print(f"  Total Overlap GO Terms Evaluated (Overlap_Count >= 2):{tot_terms:>5}")
    print("-" * 78)
    print(f"  * Overlap GO Terms WITH Overlap Genes:               {with_genes:>5} ({pct:>5.1f}%)")
    print(f"  * Overlap GO Terms WITHOUT Overlap Genes:            {stats.get('terms_without_overlap_genes', 0):>5}")
    print("-" * 78)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze whether overlap genes (in >= 2 contrasts) contribute to overlap GO terms.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--gaf-file",
        type=str,
        default="QuickGO-annotations-1764241851574-20251127.gaf",
        help="Path to GAF annotation file.",
    )
    parser.add_argument(
        "--logor",
        type=float,
        default=0.5,
        help="Log Odds Ratio threshold used in ORA overlap filenames.",
    )
    parser.add_argument(
        "--lfc",
        "--log2fc",
        dest="lfc",
        type=float,
        default=0.5,
        help="Log2 fold-change threshold used in FE overlap gene filenames.",
    )
    parser.add_argument(
        "--fe-dir",
        type=str,
        default=".",
        help="Directory containing FE overlap gene CSV files.",
    )
    parser.add_argument(
        "--ora-dir",
        type=str,
        default=".",
        help="Directory containing ORA overlap term CSV files.",
    )
    parser.add_argument(
        "--min-gene-overlap",
        type=int,
        default=2,
        help="Minimum contrast count for a gene to be considered an overlap gene.",
    )
    parser.add_argument(
        "--min-term-overlap",
        type=int,
        default=2,
        help="Minimum contrast count for a GO term to be considered an overlap term.",
    )
    parser.add_argument(
        "--all-terms",
        action="store_true",
        help="If set, evaluate all GO terms in the ORA overlap files instead of only Overlap_Count >= min_term_overlap.",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        type=str,
        default=".",
        help="Directory where contribution result CSV files will be saved.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="go_term_overlap_gene_contribution",
        help="Filename prefix for output CSV files.",
    )

    args = parser.parse_args()

    base_dir = Path.cwd()
    gaf_path = (base_dir / args.gaf_file).resolve()
    fe_dir = (base_dir / args.fe_dir).resolve()
    ora_dir = (base_dir / args.ora_dir).resolve()
    out_dir = (base_dir / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print("==============================================================================")
    print("GO TERM OVERLAP GENE CONTRIBUTION ANALYSIS")
    print("==============================================================================")
    print(f"GAF Annotation File: {gaf_path}")
    print(f"LogOR Threshold:     {args.logor}")
    print(f"Gene Log2FC Thresh:  {args.lfc}")
    print(f"Min Gene Overlap:    >= {args.min_gene_overlap} contrasts")
    print(f"Min Term Overlap:    >= {args.min_term_overlap} contrasts (all_terms={args.all_terms})")
    print(f"Output Directory:    {out_dir}")
    print("==============================================================================")

    # 1. Parse GAF
    print(f"\n[+] Loading GAF annotation file...")
    try:
        term2genes_gaf, gene2name = parse_gaf_annotations(gaf_path)
        print(f"    - Parsed {len(term2genes_gaf)} GO terms with valid gene annotations.")
    except Exception as e:
        print(f"[-] Error loading GAF file: {e}", file=sys.stderr)
        sys.exit(1)

    logor_str = format_threshold(args.logor)
    lfc_str = format_threshold(args.lfc)

    directions = [
        ("up", f"fe_overlap_up_lfc{lfc_str}.csv", f"ora_overlap_up_logor{logor_str}.csv"),
        ("down", f"fe_overlap_down_lfc{lfc_str}.csv", f"ora_overlap_down_logor{logor_str}.csv"),
        ("combined", f"fe_overlap_combined_lfc{lfc_str}.csv", f"ora_overlap_all_logor{logor_str}.csv"),
    ]

    saved_files = []

    for d, fe_name, ora_name in directions:
        fe_file = fe_dir / fe_name
        ora_file = ora_dir / ora_name

        res_df, stats = analyze_direction(
            direction=d,
            fe_file=fe_file,
            ora_file=ora_file,
            term2genes_gaf=term2genes_gaf,
            min_gene_overlap=args.min_gene_overlap,
            min_term_overlap=args.min_term_overlap,
            all_terms=args.all_terms,
        )

        print_contribution_summary(d, stats, args.lfc, args.logor)

        out_name = f"{args.prefix}_{d}_logor{logor_str}_lfc{lfc_str}.csv"
        out_file = out_dir / out_name
        saved = safe_to_csv(res_df, out_file)
        if saved:
            saved_files.append(saved)
            print(f"  [+] Saved {d.upper():<8} contribution table -> {out_file.name} ({len(res_df)} rows)")

        # Preview sample terms with overlap genes
        if not res_df.empty:
            preview_top = res_df[res_df["Has_Overlap_Gene"]].head(5)
            if not preview_top.empty:
                cols_show = ["ID", "Ontology", "Description", "N_All_Overlap_Genes", "All_Overlap_Genes"]
                print("\n    Sample Overlap GO Terms with Overlap Genes:")
                print("    " + preview_top[cols_show].to_string(index=False).replace("\n", "\n    "))
            else:
                print("    (No evaluated overlap terms had annotated overlap genes under current thresholds)")

    print("\n" + "=" * 78)
    print(f"[+] All 3 contribution files successfully generated:")
    for f in saved_files:
        print(f"    * {f.name}")
    print("=" * 78 + "\n")


if __name__ == "__main__":
    main()

