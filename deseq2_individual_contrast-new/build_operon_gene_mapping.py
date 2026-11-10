#!/usr/bin/env python3
"""
build_operon_gene_mapping.py

Phase 1 Implementation:
Parses the Staphylococcus aureus NCTC 8325 operon database (list_of_operons.txt)
and resolves all gene entities to:
  1. Genome annotation (locus tags 'SAOUHSC_...' and official gene symbols) from GCF_000013425.1.gff
  2. Study expression ortholog clusters ('cluster_id') via bidirectional matching:
     - Method 1 (Direct): matching locus tag and gene symbol against All_Genes_Locus_Tags
     - Method 2 (Ortholog): bridging NCTC 8325 protein IDs (YP_...) via MMseqs2 cluster alignments

Outputs:
  - operon_gene_mapping.csv : Per-gene annotation and mapping table (2,842 rows)
  - operon_summary.csv      : Per-operon architecture and study coverage table (1,523 rows)
"""

import os
import re
import collections
import pandas as pd

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    operon_file = os.path.join(base_dir, "list_of_operons.txt")
    gff_file = os.path.join(base_dir, "GCF_000013425.1.gff")
    mmseqs_file = os.path.join(base_dir, "mmseqs2", "cluster_alignments.tsv")
    deseq2_res_file = os.path.join(base_dir, "results_sample_info-all", "DESeq2_contrast_results.csv")
    gene_info_excel = os.path.join(base_dir, "annote", "GeneSpecificInformation_NCTC8325.xlsx")

    out_gene_mapping = os.path.join(base_dir, "operon_gene_mapping.csv")
    out_operon_summary = os.path.join(base_dir, "operon_summary.csv")

    print("=== Step 1: Parsing Genome Annotation (GCF_000013425.1.gff) ===")
    gene_info = {}
    cds_to_parent = {}
    gene_to_yp = {}
    locus_to_sym = {}
    sym_to_locus = collections.defaultdict(list)

    with open(gff_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feat_type = parts[2]
            attr_str = parts[8]
            attrs = dict(re.findall(r"([^;=]+)=([^;=]*)", attr_str))

            if feat_type in ("gene", "pseudogene"):
                gid = attrs.get("ID")
                lt = attrs.get("locus_tag", "")
                sym = attrs.get("gene", attrs.get("Name", ""))
                d = {
                    "locus_tag": lt,
                    "gene_symbol": sym,
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6]
                }
                if gid:
                    gene_info[gid] = d
                if lt:
                    gene_info[lt] = d
                    locus_to_sym[lt] = sym
                    if sym and sym != lt:
                        sym_to_locus[sym.lower()].append(lt)

            elif feat_type in ("CDS", "tRNA", "rRNA"):
                cid = attrs.get("ID")
                parent = attrs.get("Parent")
                yp = attrs.get("protein_id")
                lt = attrs.get("locus_tag")
                sym = attrs.get("gene")

                if cid and parent:
                    cds_to_parent[cid] = parent
                    if cid.startswith("cds-"):
                        cds_to_parent[cid[4:]] = parent

                if parent and yp:
                    gene_to_yp[parent] = yp
                if cid and yp:
                    gene_to_yp[cid] = yp
                    if cid.startswith("cds-"):
                        gene_to_yp[cid[4:]] = yp
                if lt and yp:
                    gene_to_yp[lt] = yp
                if sym and yp:
                    gene_to_yp[sym.lower()] = yp
                if lt and sym:
                    locus_to_sym[lt] = sym
                    sym_to_locus[sym.lower()].append(lt)

    # Supplement gene symbols and synonyms from GeneSpecificInformation_NCTC8325.xlsx if available
    if os.path.exists(gene_info_excel):
        print("  - Supplementing symbols from GeneSpecificInformation_NCTC8325.xlsx...")
        df_nctc = pd.read_excel(gene_info_excel)
        for _, row in df_nctc.iterrows():
            lt = str(row.get("locus tag", "")).strip()
            sym = str(row.get("symbol", "")).strip()
            syn = str(row.get("synonym", "")).strip()
            if lt and lt.startswith("SAOUHSC_"):
                if sym and sym != "nan" and not sym.startswith("SAOUHSC_"):
                    if lt not in locus_to_sym or not locus_to_sym[lt] or locus_to_sym[lt] == lt:
                        locus_to_sym[lt] = sym
                    sym_to_locus[sym.lower()].append(lt)
                if syn and syn != "nan":
                    sym_to_locus[syn.lower()].append(lt)

    # Specific paralog duplicates in DOOR with suffix '_2'
    paralog_fixes = {
        "rpmG_2": ("SAOUHSC_01651", "rpmG", "YP_500163.1"),
        "cbiO_2": ("SAOUHSC_02483", "cbiO", "YP_500951.1"),
        "rpsN_2": ("SAOUHSC_02499", "rpsN", "YP_500966.2"),
        "secY_2": ("SAOUHSC_02989", "secY", "YP_501438.1")
    }

    print(f"  ✓ Indexed {len(locus_to_sym)} locus tags and {len(sym_to_locus)} gene symbols from NCTC 8325.")

    print("\n=== Step 2: Indexing Study Universe (DESeq2_contrast_results.csv) ===")
    df_res = pd.read_csv(deseq2_res_file)
    study_cids = set(df_res["cluster_id"])
    print(f"  ✓ Total study clusters in meta-analysis: {len(study_cids)}")

    # Index All_Genes_Locus_Tags tokens -> set of cluster_ids
    token_to_cids = collections.defaultdict(set)
    for _, row in df_res.iterrows():
        cid = row["cluster_id"]
        raw_tags = str(row["All_Genes_Locus_Tags"])
        rep = str(row.get("Representative_Gene", ""))
        toks = [t.strip() for t in raw_tags.split("/") if t.strip()]
        if rep and rep != "nan":
            toks.extend([t.strip() for t in rep.split("/") if t.strip()])
        for t in set(toks):
            token_to_cids[t].add(cid)
            token_to_cids[t.lower()].add(cid)

    print(f"  ✓ Indexed {len(token_to_cids)} distinct identifier tokens in All_Genes_Locus_Tags.")

    print("\n=== Step 3: Indexing MMseqs2 Ortholog Alignments ===")
    mmseqs_map = {}
    if os.path.exists(mmseqs_file):
        with open(mmseqs_file, "r", encoding="utf-8") as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 2 and p[1] not in mmseqs_map:
                    mmseqs_map[p[1]] = p[0]
        print(f"  ✓ Loaded MMseqs2 ortholog mapping for {len(mmseqs_map)} proteins.")
    else:
        print("  ⚠ MMseqs2 file not found, proceeding without ortholog alignments.")

    print("\n=== Step 4: Parsing Operons and Linking to Study Clusters ===")
    operon_records = []
    with open(operon_file, "r", encoding="utf-8") as f:
        f.readline() # Header
        cur_op = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) == 1:
                try:
                    cur_op = int(parts[0])
                    continue
                except ValueError:
                    pass
            if parts[0] == "":
                parts = parts[1:]

            gid = parts[0]
            gtype = parts[1] if len(parts) > 1 else "CDS"
            cog = parts[2] if len(parts) > 2 else ""
            pos_l = int(parts[3]) if len(parts) > 3 else None
            pos_r = int(parts[4]) if len(parts) > 4 else None
            strand = parts[5] if len(parts) > 5 else ""
            func = parts[6] if len(parts) > 6 else ""

            # Resolve NCTC 8325 locus tag, symbol, and protein ID
            lt, sym, yp = None, None, None
            if gid in paralog_fixes:
                lt, sym, yp = paralog_fixes[gid]
            elif gid in cds_to_parent:
                parent = cds_to_parent[gid]
                info = gene_info.get(parent, {})
                lt = info.get("locus_tag")
                sym = locus_to_sym.get(lt, info.get("gene_symbol"))
                yp = gene_to_yp.get(gid, gene_to_yp.get(parent))
            elif gid in gene_info:
                info = gene_info[gid]
                lt = info.get("locus_tag")
                sym = locus_to_sym.get(lt, info.get("gene_symbol"))
                yp = gene_to_yp.get(gid, gene_to_yp.get(lt))
            elif gid.lower() in sym_to_locus:
                lt = sym_to_locus[gid.lower()][0]
                sym = locus_to_sym.get(lt, gid)
                yp = gene_to_yp.get(gid.lower(), gene_to_yp.get(lt))

            # Clean up symbol if it's identical to locus_tag
            clean_sym = sym if (sym and sym != lt) else None

            # ----------------------------------------------------
            # Method 1 (Direct User Request):
            # Match using locus tag and gene symbol in All_Genes_Locus_Tags
            # ----------------------------------------------------
            m1_cids = set()
            m1_sources = []

            # Match locus tag
            if lt and lt in token_to_cids:
                m1_cids.update(token_to_cids[lt])
                m1_sources.append("locus_tag")

            # Match gene symbol
            if clean_sym:
                if clean_sym in token_to_cids:
                    m1_cids.update(token_to_cids[clean_sym])
                    m1_sources.append("gene_symbol")
                elif clean_sym.lower() in token_to_cids:
                    m1_cids.update(token_to_cids[clean_sym.lower()])
                    m1_sources.append("gene_symbol")

            # Match original operon gene id if it represents a symbol
            if gid and gid not in (lt, clean_sym):
                if gid in token_to_cids or gid.lower() in token_to_cids:
                    c = token_to_cids.get(gid, token_to_cids.get(gid.lower(), set()))
                    m1_cids.update(c)
                    m1_sources.append("operon_symbol")

            direct_cid = list(m1_cids)[0] if m1_cids else None
            direct_source = "/".join(m1_sources) if m1_sources else "unmapped"

            # ----------------------------------------------------
            # Method 2 (Ortholog MMseqs2):
            # Bridge YP_ protein to cluster_id
            # ----------------------------------------------------
            ortholog_cid = None
            if yp and yp in mmseqs_map:
                c = mmseqs_map[yp]
                if c in study_cids:
                    ortholog_cid = c

            # ----------------------------------------------------
            # Combined Consensus Resolution:
            # Prefer direct locus/symbol match; supplement with ortholog cluster
            # ----------------------------------------------------
            if direct_cid:
                final_cid = direct_cid
                final_source = direct_source
            elif ortholog_cid:
                final_cid = ortholog_cid
                final_source = "ortholog_cluster"
            else:
                final_cid = None
                final_source = "unmapped"

            in_study = (final_cid in study_cids) if final_cid else False

            operon_records.append({
                "operon_id": cur_op,
                "operon_gene_id": gid,
                "gene_type": gtype,
                "locus_tag": lt or "",
                "gene_symbol": clean_sym or "",
                "protein_id": yp or "",
                "pos_left": pos_l,
                "pos_right": pos_r,
                "strand": strand,
                "cog_id": cog,
                "function": func,
                "direct_cluster_id": direct_cid or "",
                "direct_match_source": direct_source,
                "ortholog_cluster_id": ortholog_cid or "",
                "final_cluster_id": final_cid or "",
                "final_match_source": final_source,
                "in_study_universe": in_study
            })

    df_gene_map = pd.DataFrame(operon_records)
    df_gene_map.to_csv(out_gene_mapping, index=False)
    print(f"  ✓ Exported detailed gene mapping to: {out_gene_mapping}")

    # ----------------------------------------------------
    # Step 5: Build Operon-Level Summary Table
    # ----------------------------------------------------
    print("\n=== Step 5: Generating Operon Architecture Summary ===")
    operon_summaries = []

    for op_id, grp in df_gene_map.groupby("operon_id"):
        tot_genes = len(grp)
        mapped_grp = grp[grp["in_study_universe"] & (grp["final_cluster_id"] != "")]
        mapped_genes_count = len(mapped_grp)
        mapped_cids = sorted(list(set(mapped_grp["final_cluster_id"])))

        symbols = [s for s in grp["gene_symbol"] if s]
        locus_tags = [lt for lt in grp["locus_tag"] if lt]

        # Operon representative label
        if symbols:
            op_label = f"Operon_{op_id} ({'/'.join(symbols[:4])}{'...' if len(symbols)>4 else ''})"
        elif locus_tags:
            op_label = f"Operon_{op_id} ({'/'.join(locus_tags[:3])}{'...' if len(locus_tags)>3 else ''})"
        else:
            op_label = f"Operon_{op_id}"

        # COG functional classification
        cogs = [c for c in grp["cog_id"] if c and c != "NA"]
        dom_cog = cogs[0] if cogs else "Unassigned"

        funcs = [f for f in grp["function"] if f and f != "NA"]
        prim_func = funcs[0] if funcs else "Uncharacterized"

        strand = grp["strand"].iloc[0]
        pos_min = grp["pos_left"].min()
        pos_max = grp["pos_right"].max()

        operon_summaries.append({
            "operon_id": op_id,
            "operon_label": op_label,
            "strand": strand,
            "pos_left": pos_min,
            "pos_right": pos_max,
            "span_bp": pos_max - pos_min + 1 if (pos_min and pos_max) else 0,
            "total_genes": tot_genes,
            "is_polycistronic": tot_genes >= 2,
            "mapped_study_genes": mapped_genes_count,
            "coverage_ratio": round(mapped_genes_count / tot_genes, 4),
            "study_clusters_count": len(mapped_cids),
            "study_cluster_ids": "/".join(mapped_cids),
            "gene_symbols": "/".join(symbols),
            "locus_tags": "/".join(locus_tags),
            "dominant_cog": dom_cog,
            "primary_function": prim_func
        })

    df_op_sum = pd.DataFrame(operon_summaries)
    df_op_sum.to_csv(out_operon_summary, index=False)
    print(f"  ✓ Exported operon summary to: {out_operon_summary}")

    # Summary Statistics
    print("\n=== Phase 1 Completion Statistics ===")
    print(f"Total Operons: {len(df_op_sum)}")
    print(f"  ├── Polycistronic (>= 2 genes): {(df_op_sum['is_polycistronic']).sum()}")
    print(f"  └── Monocistronic (1 gene):     {(~df_op_sum['is_polycistronic']).sum()}")
    print(f"Total Genes in Operon Database: {len(df_gene_map)}")
    print(f"  ├── Mapped directly via All_Genes_Locus_Tags (locus+symbol): {(df_gene_map['direct_cluster_id'] != '').sum()} ({(df_gene_map['direct_cluster_id'] != '').sum()/len(df_gene_map)*100:.1f}%)")
    print(f"  └── Mapped via Full Consensus (Direct + Ortholog):           {(df_gene_map['in_study_universe']).sum()} ({(df_gene_map['in_study_universe']).sum()/len(df_gene_map)*100:.1f}%)")
    print(f"Study Clusters Covered in DESeq2 Meta-Analysis: {df_gene_map['final_cluster_id'].replace('', pd.NA).dropna().nunique()} / {len(study_cids)} ({df_gene_map['final_cluster_id'].replace('', pd.NA).dropna().nunique()/len(study_cids)*100:.1f}%)")
    print(f"Operons with >= 1 gene in study universe: {(df_op_sum['mapped_study_genes'] >= 1).sum()} / {len(df_op_sum)} ({(df_op_sum['mapped_study_genes'] >= 1).sum()/len(df_op_sum)*100:.1f}%)")
    print(f"Operons with >= 50% genes in study universe: {(df_op_sum['coverage_ratio'] >= 0.5).sum()} / {len(df_op_sum)} ({(df_op_sum['coverage_ratio'] >= 0.5).sum()/len(df_op_sum)*100:.1f}%)")

if __name__ == "__main__":
    main()

