#!/usr/bin/env python3
"""
build_operon_gene_mapping.py

Phase 1 Implementation:
Parses the Staphylococcus aureus NCTC 8325 operon database (list_of_operons.txt)
and resolves all gene entities to:
  1. Genome annotation (locus tags 'SAOUHSC_...' and official gene symbols) from GCF_000013425.1.gff
  2. AureoWiki multi-strain pan-genome annotations (pan gene symbols, synonyms, and SAUPAN ortholog tokens)
  3. Study expression ortholog clusters ('cluster_id') via multi-method matching:
     - Method 1 (Direct): matching locus tag and gene symbol against All_Genes_Locus_Tags
     - Method 2 (AureoWiki Ortholog): bridging via AureoWiki SAUPAN pan-locus groups and pan gene symbols
     - Method 3 (MMseqs2 Ortholog): bridging NCTC 8325 protein IDs (YP_...) via MMseqs2 cluster alignments

Outputs:
  - operon_gene_mapping.csv : Per-gene annotation and mapping table (2,842 rows)
  - operon_summary.csv      : Per-operon architecture and study coverage table (1,523 rows)
"""

import os
import re
import glob
import collections
import pandas as pd
import networkx as nx

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    operon_file = os.path.join(base_dir, "list_of_operons.txt")
    gff_file = os.path.join(base_dir, "GCF_000013425.1.gff")
    mmseqs_file = os.path.join(base_dir, "mmseqs2", "cluster_alignments.tsv")
    annote_dir = os.path.join(base_dir, "annote")

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

    # Step 1b: Supplement symbols and cross-strain orthologs from AureoWiki tables
    strain_files = sorted(glob.glob(os.path.join(annote_dir, "GeneSpecificInformation_*.xlsx")))
    pan_lt_to_tokens = collections.defaultdict(set)
    saouhsc_to_pan_lt = {}

    if strain_files:
        print(f"  - Supplementing annotations from {len(strain_files)} AureoWiki strain files...")
        for sf in strain_files:
            try:
                df_strain = pd.read_excel(sf)
            except Exception as e:
                print(f"    ⚠ Could not read {sf}: {e}")
                continue

            for _, row in df_strain.iterrows():
                lt = str(row.get("locus tag", "")).strip()
                pan_lt = str(row.get("pan locus tag", "")).strip()
                pan_sym = str(row.get("pan gene symbol", "")).strip()
                sym = str(row.get("symbol", "")).strip()
                syn = str(row.get("synonym", "")).strip()
                valid_pan_lt = pan_lt if (pan_lt and pan_lt != "nan" and pan_lt.startswith("SAUPAN")) else None

                cur_sym = None
                if pan_sym and pan_sym != "nan" and not pan_sym.startswith("SAOUHSC_") and not pan_sym.startswith("SA_") and not pan_sym.startswith("SAUSA300_"):
                    cur_sym = pan_sym
                elif sym and sym != "nan" and not sym.startswith("SAOUHSC_") and not sym.startswith("SA_") and not sym.startswith("SAUSA300_"):
                    cur_sym = sym

                if lt.startswith("SAOUHSC_"):
                    if cur_sym and (lt not in locus_to_sym or not locus_to_sym[lt] or locus_to_sym[lt] == lt):
                        locus_to_sym[lt] = cur_sym
                        sym_to_locus[cur_sym.lower()].append(lt)
                    if valid_pan_lt:
                        saouhsc_to_pan_lt[lt] = valid_pan_lt

                if valid_pan_lt:
                    if lt and lt != "nan":
                        pan_lt_to_tokens[valid_pan_lt].add(lt)
                    if cur_sym:
                        pan_lt_to_tokens[valid_pan_lt].add(cur_sym)
                    if syn and syn != "nan":
                        for s in str(syn).split(";"):
                            s = s.strip()
                            if s and not s.startswith("SAOUHSC_") and not s.startswith("SA_") and not s.startswith("SAUSA300_"):
                                pan_lt_to_tokens[valid_pan_lt].add(s)

    # Specific paralog duplicates in DOOR with suffix '_2'
    paralog_fixes = {
        "rpmG_2": ("SAOUHSC_01651", "rpmG", "YP_500163.1"),
        "cbiO_2": ("SAOUHSC_02483", "cbiO", "YP_500951.1"),
        "rpsN_2": ("SAOUHSC_02499", "rpsN", "YP_500966.2"),
        "secY_2": ("SAOUHSC_02989", "secY", "YP_501438.1")
    }

    print(f"  ✓ Indexed {len(locus_to_sym)} locus tags and {len(sym_to_locus)} gene symbols from NCTC 8325 + AureoWiki.")
    print(f"  ✓ Indexed {len(pan_lt_to_tokens)} AureoWiki pan-locus (SAUPAN) orthology groups.")

    print("\n=== Step 2: Indexing Study Universe (DESeq2_contrast_results.csv) ===")
    deseq2_files = sorted(glob.glob(os.path.join(base_dir, "results_*", "DESeq2_contrast_results.csv")))
    if not deseq2_files:
        deseq2_files = [os.path.join(base_dir, "results_sample_info-all", "DESeq2_contrast_results.csv")]

    study_cids = set()
    token_to_cids = collections.defaultdict(set)

    for d_file in deseq2_files:
        if not os.path.exists(d_file):
            continue
        df_res = pd.read_csv(d_file)
        study_cids.update(df_res["cluster_id"].dropna().astype(str).str.strip())
        for _, row in df_res.iterrows():
            cid = str(row["cluster_id"]).strip()
            raw_tags = str(row.get("All_Genes_Locus_Tags", ""))
            rep = str(row.get("Representative_Gene", ""))
            toks = [t.strip() for t in raw_tags.split("/") if t.strip()]
            if rep and rep != "nan":
                toks.extend([t.strip() for t in rep.split("/") if t.strip()])
            for t in set(toks):
                token_to_cids[t].add(cid)
                token_to_cids[t.lower()].add(cid)

    print(f"  ✓ Total unique study clusters across results: {len(study_cids)}")
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
            # Method 1 (Direct Locus & Symbol Matching):
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

            # Match cross-strain ortholog tokens via AureoWiki SAUPAN group
            if lt and lt in saouhsc_to_pan_lt:
                pan_id = saouhsc_to_pan_lt[lt]
                for orth_tok in pan_lt_to_tokens[pan_id]:
                    if orth_tok in token_to_cids or orth_tok.lower() in token_to_cids:
                        c = token_to_cids.get(orth_tok, token_to_cids.get(orth_tok.lower(), set()))
                        m1_cids.update(c)
                        if "aureowiki_ortholog" not in m1_sources:
                            m1_sources.append("aureowiki_ortholog")

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
            # Prefer direct locus/symbol/AureoWiki match; supplement with MMseqs2
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

    # ----------------------------------------------------
    # Step 4b: Alternative 4 - Operon Family Merging
    # Whenever >= 2 DOOR operons share an expression cluster (paralogs, duplicates),
    # merge them into a single non-redundant operon family unit.
    # ----------------------------------------------------
    print("\n=== Step 4b: Merging Multi-Operon Clusters into Non-Redundant Operon Families (Alt 4) ===")
    G = nx.Graph()
    for oid in df_gene_map["operon_id"].unique():
        G.add_node(oid)

    valid_cids = df_gene_map[df_gene_map["final_cluster_id"].notna() & (df_gene_map["final_cluster_id"] != "")]
    for cid, grp in valid_cids.groupby("final_cluster_id"):
        oids = grp["operon_id"].unique()
        if len(oids) > 1:
            for i in range(len(oids)):
                for j in range(i + 1, len(oids)):
                    G.add_edge(oids[i], oids[j])

    components = sorted(list(nx.connected_components(G)), key=lambda c: min(c))
    multi_count = sum(1 for c in components if len(c) > 1)
    print(f"  Identified {multi_count} multi-operon family components across {len(G.nodes)} original operons.")
    print(f"  Total non-redundant operon units after merging: {len(components)}")

    comp_map = {}
    is_merged_map = {}
    orig_ops_map = {}

    for comp in components:
        sorted_oids = sorted(list(comp))
        if len(sorted_oids) == 1:
            m_id = str(sorted_oids[0])
            is_m = False
            o_str = str(sorted_oids[0])
        else:
            m_id = "_".join(str(o) for o in sorted_oids)
            is_m = True
            o_str = "/".join(str(o) for o in sorted_oids)
        for o in sorted_oids:
            comp_map[o] = m_id
            is_merged_map[o] = is_m
            orig_ops_map[o] = o_str

    df_gene_map["original_operon_id"] = df_gene_map["operon_id"]
    df_gene_map["operon_id"] = df_gene_map["original_operon_id"].map(comp_map)
    df_gene_map["is_merged_operon"] = df_gene_map["original_operon_id"].map(is_merged_map)
    df_gene_map["original_operon_ids"] = df_gene_map["original_operon_id"].map(orig_ops_map)

    # Reorder columns
    cols = [
        "operon_id", "original_operon_id", "is_merged_operon", "original_operon_ids",
        "operon_gene_id", "gene_type", "locus_tag", "gene_symbol", "protein_id",
        "pos_left", "pos_right", "strand", "cog_id", "function",
        "direct_cluster_id", "direct_match_source", "ortholog_cluster_id",
        "final_cluster_id", "final_match_source", "in_study_universe"
    ]
    df_gene_map = df_gene_map[cols]
    df_gene_map.to_csv(out_gene_mapping, index=False)
    print(f"  ✓ Exported detailed non-redundant gene mapping to: {out_gene_mapping}")

    # ----------------------------------------------------
    # Step 5: Build Operon-Level Summary Table
    # ----------------------------------------------------
    print("\n=== Step 5: Generating Operon Architecture Summary ===")
    operon_summaries = []

    for op_id, grp in df_gene_map.groupby("operon_id", sort=False):
        tot_genes = len(grp)
        mapped_grp = grp[grp["in_study_universe"] & (grp["final_cluster_id"] != "")]
        mapped_genes_count = len(mapped_grp)
        mapped_cids = sorted(list(set(mapped_grp["final_cluster_id"])))

        symbols = [s for s in grp["gene_symbol"] if s and s != "-"]
        unique_symbols = list(dict.fromkeys(symbols))

        locus_tags = [lt for lt in grp["locus_tag"] if lt]
        unique_locus_tags = list(dict.fromkeys(locus_tags))

        is_merged = grp["is_merged_operon"].iloc[0]
        orig_ops = grp["original_operon_ids"].iloc[0]

        # Operon representative label
        if unique_symbols:
            op_label = f"Operon_{op_id} ({'/'.join(unique_symbols[:4])}{'...' if len(unique_symbols)>4 else ''})"
        elif unique_locus_tags:
            op_label = f"Operon_{op_id} ({'/'.join(unique_locus_tags[:3])}{'...' if len(unique_locus_tags)>3 else ''})"
        else:
            op_label = f"Operon_{op_id}"

        # COG functional classification
        cogs = [c for c in grp["cog_id"] if c and c != "NA"]
        dom_cog = cogs[0] if cogs else "Unassigned"

        funcs = [f for f in grp["function"] if f and f != "NA"]
        prim_func = funcs[0] if funcs else "Uncharacterized"

        strands = grp["strand"].dropna().unique()
        strand = strands[0] if len(strands) == 1 else "mixed"
        pos_min = grp["pos_left"].min()
        pos_max = grp["pos_right"].max()

        operon_summaries.append({
            "operon_id": op_id,
            "original_operon_ids": orig_ops,
            "is_merged_operon": is_merged,
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
            "gene_symbols": "/".join(unique_symbols),
            "locus_tags": "/".join(unique_locus_tags),
            "dominant_cog": dom_cog,
            "primary_function": prim_func
        })

    operon_summaries.sort(key=lambda x: int(str(x["operon_id"]).split("_")[0]))
    df_op_sum = pd.DataFrame(operon_summaries)
    df_op_sum.to_csv(out_operon_summary, index=False)
    print(f"  ✓ Exported non-redundant operon summary to: {out_operon_summary}")

    # Summary Statistics
    print("\n=== Phase 1 Completion Statistics ===")
    print(f"Original DOOR Operons: {len(G.nodes)}")
    print(f"Non-Redundant Operon Units: {len(df_op_sum)}")
    print(f"  ├── Single Unmerged Operons: {(~df_op_sum['is_merged_operon']).sum()}")
    print(f"  └── Multi-Operon Merged Families: {(df_op_sum['is_merged_operon']).sum()}")
    print(f"Total Genes in Operon Database: {len(df_gene_map)}")
    print(f"  ├── Mapped directly via All_Genes_Locus_Tags (locus+symbol+AureoWiki): {(df_gene_map['direct_cluster_id'] != '').sum()} ({(df_gene_map['direct_cluster_id'] != '').sum()/len(df_gene_map)*100:.1f}%)")
    print(f"  └── Mapped via Full Consensus (Direct + Ortholog):                     {(df_gene_map['in_study_universe']).sum()} ({(df_gene_map['in_study_universe']).sum()/len(df_gene_map)*100:.1f}%)")
    print(f"Study Clusters Covered in DESeq2 Meta-Analysis: {df_gene_map['final_cluster_id'].replace('', pd.NA).dropna().nunique()} / {len(study_cids)} ({df_gene_map['final_cluster_id'].replace('', pd.NA).dropna().nunique()/len(study_cids)*100:.1f}%)")
    print(f"Operons with >= 1 gene in study universe: {(df_op_sum['mapped_study_genes'] >= 1).sum()} / {len(df_op_sum)} ({(df_op_sum['mapped_study_genes'] >= 1).sum()/len(df_op_sum)*100:.1f}%)")
    print(f"Operons with >= 50% genes in study universe: {(df_op_sum['coverage_ratio'] >= 0.5).sum()} / {len(df_op_sum)} ({(df_op_sum['coverage_ratio'] >= 0.5).sum()/len(df_op_sum)*100:.1f}%)")

if __name__ == "__main__":
    main()
