# ==============================================================================
# Script: gene_id_mapping.R
# Purpose: Standardized cluster_id <-> locus-tag/gene mapping for iModulon
#          analysis, mimicking the 3-method standardization pipeline in
#          build_operon_gene_mapping.py (without altering ora_meta_analysis.R):
#            Method 1 (Direct): GCF_000013425.1.gff + iModulonDB gene_info/trn
#            Method 2 (AureoWiki Ortholog): SAUPAN groups, symbols & synonyms
#            Method 3 (MMseqs2 & GAF/CI fallback): protein_id -> cluster_id
# ==============================================================================

#' Count leading GAF header lines starting with '!'
#' @keywords internal
.count_gaf_header_lines <- function(gaf_file) {
  con <- file(gaf_file, open = "r")
  on.exit(close(con), add = TRUE)
  n_skip <- 0L
  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (length(ln) == 0L || !startsWith(ln, "!")) break
    n_skip <- n_skip + 1L
  }
  n_skip
}

#' Expand delimited values in parallel with an ID vector (vectorized)
#' @keywords internal
.expand_token_pairs <- function(ids, raw_vals, sep_regex = "[/;,]") {
  ids <- as.character(ids)
  vals <- trimws(as.character(raw_vals))
  keep <- !is.na(ids) & ids != "" & !is.na(vals) & !vals %in% c("", "nan", "NaN", "NA", "-")
  if (!any(keep)) return(data.frame(id = character(0), token = character(0), stringsAsFactors = FALSE))

  ids <- ids[keep]
  vals <- vals[keep]
  has_sep <- grepl(sep_regex, vals)

  parts <- list()
  if (any(!has_sep)) {
    parts[[1]] <- data.frame(id = ids[!has_sep], token = vals[!has_sep], stringsAsFactors = FALSE)
  }
  if (any(has_sep)) {
    sp <- strsplit(vals[has_sep], sep_regex)
    lens <- lengths(sp)
    toks <- trimws(unlist(sp, use.names = FALSE))
    rep_ids <- rep(ids[has_sep], lens)
    ok <- !is.na(toks) & !toks %in% c("", "nan", "NaN", "NA", "-")
    if (any(ok)) {
      parts[[2]] <- data.frame(id = rep_ids[ok], token = toks[ok], stringsAsFactors = FALSE)
    }
  }

  out <- do.call(rbind, parts)
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame(id = character(0), token = character(0), stringsAsFactors = FALSE))
  }
  unique(out)
}

#' Parse NCTC 8325 GFF (GCF_000013425.1.gff) for locus_tag <-> gene_symbol <-> protein_id
#' Mirrors Step 1 of build_operon_gene_mapping.py.
#' @keywords internal
.parse_nctc8325_gff <- function(gff_file) {
  empty <- data.frame(locus_tag = character(0), gene_symbol = character(0), protein_id = character(0), stringsAsFactors = FALSE)
  if (is.null(gff_file) || !file.exists(gff_file)) return(empty)

  lines <- readLines(gff_file, warn = FALSE)
  lines <- lines[!startsWith(lines, "#") & nzchar(lines)]
  if (length(lines) == 0) return(empty)

  parts <- strsplit(lines, "\t", fixed = TRUE)
  valid <- lengths(parts) >= 9L
  parts <- parts[valid]
  if (length(parts) == 0) return(empty)

  feat_types <- vapply(parts, `[`, character(1), 3L)
  keep <- feat_types %in% c("gene", "pseudogene", "CDS", "tRNA", "rRNA")
  if (!any(keep)) return(empty)

  attrs <- vapply(parts[keep], `[`, character(1), 9L)
  extract_attr <- function(key) {
    pat <- paste0("(^|;)", key, "=([^;]*)")
    m <- regexec(pat, attrs, perl = TRUE)
    res <- regmatches(attrs, m)
    vapply(res, function(x) if (length(x) >= 3L) x[3L] else NA_character_, character(1))
  }

  lt <- extract_attr("locus_tag")
  sym <- extract_attr("gene")
  yp <- extract_attr("protein_id")
  ok <- !is.na(lt) & lt != ""
  if (!any(ok)) return(empty)

  unique(data.frame(
    locus_tag = lt[ok],
    gene_symbol = ifelse(!is.na(sym[ok]) & sym[ok] != "" & sym[ok] != lt[ok], sym[ok], NA_character_),
    protein_id = ifelse(!is.na(yp[ok]) & yp[ok] != "", yp[ok], NA_character_),
    stringsAsFactors = FALSE
  ))
}

#' Build tiered pipeline cluster_id <-> iModulon target (USA300HOU_RS*) dictionary
#'
#' Mimics build_operon_gene_mapping.py's 3-method standardization pipeline:
#'   1. Direct locus/symbol matching (All_Genes_Locus_Tags, GFF, gene_info, trn)
#'   2. AureoWiki pan-genome SAUPAN orthologs (symbols, synonyms, locus tags, OrthologueTable)
#'   3. MMseqs2 cluster_alignments.tsv protein_id bridge + GAF/case-insensitive fallback
#'
#' @param cluster_df Data frame containing cluster_id, Representative_Gene, All_Genes_Locus_Tags.
#' @param m_gene_ids Character vector of row names from iModulonDB M matrix.
#' @param gene_info_file Path to staph_precise165-gene_info.csv.
#' @param trn_file Path to staph_precise165-trn.csv.
#' @param gff_file Path to NCTC 8325 reference GFF (GCF_000013425.1.gff).
#' @param mmseqs_file Path to mmseqs2/cluster_alignments.tsv.
#' @param orthologue_file Path to OrthologueTable.xlsx.
#' @param annote_dir Directory containing GeneSpecificInformation_*.xlsx files.
#' @param operon_file Path to operon_gene_mapping.csv.
#' @param gaf_file Path to QuickGO GAF file for UniProt-to-symbol bridging.
#' @param verbose Logical; print detailed coverage breakdown.
#' @return List with mapping_df (cluster_id, target, match_tier), cluster_summary, and stats.
map_clusters_to_imodulon <- function(cluster_df,
                                     m_gene_ids,
                                     gene_info_file = "imodulon_data/staph_precise165-gene_info.csv",
                                     trn_file = "imodulon_data/staph_precise165-trn.csv",
                                     gff_file = "GCF_000013425.1.gff",
                                     mmseqs_file = file.path("mmseqs2", "cluster_alignments.tsv"),
                                     orthologue_file = "OrthologueTable.xlsx",
                                     annote_dir = "annote",
                                     operon_file = "operon_gene_mapping.csv",
                                     gaf_file = "QuickGO-annotations-1764241851574-20251127.gaf",
                                     verbose = TRUE) {
  m_gene_ids <- unique(as.character(m_gene_ids))
  cluster_df <- cluster_df[!duplicated(cluster_df$cluster_id), , drop = FALSE]
  cluster_ids <- as.character(cluster_df$cluster_id)

  # --- 1. Index study cluster tokens (Step 2 of build_operon_gene_mapping.py) ---
  c_parts <- list(
    data.frame(id = cluster_ids, token = cluster_ids, stringsAsFactors = FALSE),
    data.frame(id = cluster_ids, token = sub("\\.[0-9]+$", "", cluster_ids), stringsAsFactors = FALSE)
  )
  if ("Representative_Gene" %in% colnames(cluster_df)) {
    c_parts[[length(c_parts) + 1]] <- .expand_token_pairs(cluster_ids, cluster_df$Representative_Gene, "/")
  }
  if ("All_Genes_Locus_Tags" %in% colnames(cluster_df)) {
    c_parts[[length(c_parts) + 1]] <- .expand_token_pairs(cluster_ids, cluster_df$All_Genes_Locus_Tags, "/")
  }
  if (!is.null(operon_file) && file.exists(operon_file)) {
    op_df <- tryCatch(read.csv(operon_file, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(op_df) && "final_cluster_id" %in% colnames(op_df)) {
      for (col in intersect(c("locus_tag", "gene_symbol", "protein_id"), colnames(op_df))) {
        c_parts[[length(c_parts) + 1]] <- .expand_token_pairs(op_df$final_cluster_id, op_df[[col]], "[/;,]")
      }
    }
  }
  cluster_tok_df <- unique(do.call(rbind, c_parts))
  cluster_tok_df <- cluster_tok_df[cluster_tok_df$id %in% cluster_ids, , drop = FALSE]
  colnames(cluster_tok_df) <- c("cluster_id", "token")
  cluster_tok_df$token_low <- tolower(cluster_tok_df$token)
  cluster_tok_df$pair_key <- paste0(cluster_tok_df$cluster_id, "\t", cluster_tok_df$token)

  # --- 2. Method 1 (Direct): M_ids + gene_info + trn + GFF ---
  t1_parts <- list(data.frame(id = m_gene_ids, token = m_gene_ids, stringsAsFactors = FALSE))
  gi_sub <- NULL
  if (!is.null(gene_info_file) && file.exists(gene_info_file)) {
    gi <- read.csv(gene_info_file, stringsAsFactors = FALSE)
    gi_sub <- gi[gi$gene_id %in% m_gene_ids, , drop = FALSE]
    if (nrow(gi_sub) > 0) {
      gi_sub$saouhsc <- ifelse(!is.na(gi_sub$string_id) & gi_sub$string_id != "",
                               sub("^93061\\.", "", gi_sub$string_id), NA_character_)
      for (col in intersect(c("gene_name", "saouhsc", "uniprot_id", "gene_synonyms"), colnames(gi_sub))) {
        t1_parts[[length(t1_parts) + 1]] <- .expand_token_pairs(gi_sub$gene_id, gi_sub[[col]], "[/;,|]")
      }
    }
  }

  gff_df <- .parse_nctc8325_gff(gff_file)
  if (!is.null(gi_sub) && nrow(gff_df) > 0) {
    gi_gff <- unique(merge(
      gi_sub[!is.na(gi_sub$saouhsc) & gi_sub$saouhsc != "", c("gene_id", "saouhsc")],
      gff_df,
      by.x = "saouhsc",
      by.y = "locus_tag"
    ))
    if (nrow(gi_gff) > 0) {
      t1_parts[[length(t1_parts) + 1]] <- .expand_token_pairs(gi_gff$gene_id, gi_gff$gene_symbol, "[/;,]")
      t1_parts[[length(t1_parts) + 1]] <- .expand_token_pairs(gi_gff$gene_id, gi_gff$protein_id, "[/;,]")
    }
  }

  trn_sub <- NULL
  if (!is.null(trn_file) && file.exists(trn_file)) {
    trn <- read.csv(trn_file, stringsAsFactors = FALSE)
    if (all(c("target_id", "N315_gene_id") %in% colnames(trn))) {
      trn_sub <- trn[trn$target_id %in% m_gene_ids & !is.na(trn$N315_gene_id) & trn$N315_gene_id != "", , drop = FALSE]
      if (nrow(trn_sub) > 0) {
        # Exclude 1-to-many corrupted N315_gene_id entries (e.g. SA2480 mapped to 177 targets in TRN)
        trn_pairs <- unique(trn_sub[, c("target_id", "N315_gene_id"), drop = FALSE])
        n315_counts <- table(trn_pairs$N315_gene_id)
        valid_n315 <- names(n315_counts)[n315_counts == 1L]
        trn_sub <- trn_sub[trn_sub$N315_gene_id %in% valid_n315, , drop = FALSE]
        if (nrow(trn_sub) > 0) {
          t1_parts[[length(t1_parts) + 1]] <- .expand_token_pairs(trn_sub$target_id, trn_sub$N315_gene_id, "[/;,|]")
        }
      }
    }
  }

  m_tok_df <- unique(do.call(rbind, t1_parts))
  colnames(m_tok_df) <- c("target", "token")
  m_tok_df$token_low <- tolower(m_tok_df$token)

  # QuickGO UniProt -> gene_symbol pairs (restricted to relevant UniProt IDs)
  gaf_uni_df <- data.frame(uniprot_id = character(0), symbol = character(0), stringsAsFactors = FALSE)
  if (!is.null(gaf_file) && file.exists(gaf_file) && !is.null(gi_sub)) {
    rel_unis <- unique(gi_sub$uniprot_id[!is.na(gi_sub$uniprot_id) & gi_sub$uniprot_id != ""])
    if (length(rel_unis) > 0) {
      gaf_raw <- tryCatch({
        if (requireNamespace("data.table", quietly = TRUE)) {
          n_skip <- .count_gaf_header_lines(gaf_file)
          data.table::fread(gaf_file, sep = "\t", header = FALSE, skip = n_skip,
                            select = c(2L, 3L), data.table = FALSE, quote = "")
        } else {
          read.delim(gaf_file, comment.char = "!", header = FALSE, stringsAsFactors = FALSE, quote = "")[, c(2L, 3L)]
        }
      }, error = function(e) NULL)
      if (!is.null(gaf_raw) && ncol(gaf_raw) >= 2) {
        gaf_uni_df <- unique(gaf_raw[gaf_raw[[1]] %in% rel_unis & !is.na(gaf_raw[[2]]) & gaf_raw[[2]] != "", , drop = FALSE])
        colnames(gaf_uni_df) <- c("uniprot_id", "symbol")
      }
    }
  }

  # --- 3. Method 2 (AureoWiki Ortholog): GeneSpecificInformation + OrthologueTable ---
  pan_parts <- list()
  if (!is.null(annote_dir) && dir.exists(annote_dir) && requireNamespace("openxlsx", quietly = TRUE)) {
    gsi_files <- list.files(annote_dir, pattern = "^GeneSpecificInformation_.*\\.xlsx$", full.names = TRUE)
    for (gf in gsi_files) {
      df <- tryCatch(openxlsx::read.xlsx(gf), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) next
      colnames(df) <- tolower(gsub("[. ]+", "_", colnames(df)))
      if (!"pan_locus_tag" %in% colnames(df)) next
      for (col in intersect(c("locus_tag", "old_locus_tag", "pan_gene_symbol", "symbol", "synonym", "protein_refseq", "uniprot"), colnames(df))) {
        pan_parts[[length(pan_parts) + 1]] <- .expand_token_pairs(df$pan_locus_tag, df[[col]], "[,;]")
      }
    }
  }
  if (!is.null(orthologue_file) && file.exists(orthologue_file) && requireNamespace("openxlsx", quietly = TRUE)) {
    ot <- tryCatch(openxlsx::read.xlsx(orthologue_file), error = function(e) NULL)
    if (!is.null(ot) && nrow(ot) > 0) {
      colnames(ot) <- trimws(colnames(ot))
      pan_col <- colnames(ot)[grepl("^pan[ .]?ID$", colnames(ot), ignore.case = TRUE)][1]
      if (!is.na(pan_col)) {
        for (col in setdiff(colnames(ot), pan_col)) {
          pan_parts[[length(pan_parts) + 1]] <- .expand_token_pairs(ot[[pan_col]], ot[[col]], "[,;]")
        }
      }
    }
  }

  pan_tok_df <- unique(do.call(rbind, pan_parts))
  if (is.null(pan_tok_df) || nrow(pan_tok_df) == 0) {
    pan_tok_df <- data.frame(pan_id = character(0), token = character(0), token_low = character(0), stringsAsFactors = FALSE)
  } else {
    colnames(pan_tok_df) <- c("pan_id", "token")
    pan_tok_df$token_low <- tolower(pan_tok_df$token)
  }

  m_tok_for_pan <- m_tok_df[, c("target", "token", "token_low")]
  if (nrow(gaf_uni_df) > 0 && !is.null(gi_sub)) {
    gi_u <- unique(merge(gi_sub[, c("gene_id", "uniprot_id")], gaf_uni_df, by = "uniprot_id"))
    if (nrow(gi_u) > 0) {
      m_tok_for_pan <- unique(rbind(
        m_tok_for_pan,
        data.frame(target = gi_u$gene_id, token = gi_u$symbol, token_low = tolower(gi_u$symbol), stringsAsFactors = FALSE)
      ))
    }
  }
  pan_to_m_df <- unique(merge(
    pan_tok_df[, c("pan_id", "token")],
    m_tok_for_pan[, c("target", "token")],
    by = "token"
  )[, c("pan_id", "target")])

  # --- 4. Method 3 (MMseqs2 Ortholog): protein_id -> cluster_id ---
  mm_cluster_df <- data.frame(cluster_id = character(0), protein_id = character(0), stringsAsFactors = FALSE)
  if (!is.null(mmseqs_file) && file.exists(mmseqs_file)) {
    mm_raw <- tryCatch({
      if (requireNamespace("data.table", quietly = TRUE)) {
        data.table::fread(mmseqs_file, sep = "\t", header = FALSE, select = c(1L, 2L), data.table = FALSE)
      } else {
        read.delim(mmseqs_file, header = FALSE, stringsAsFactors = FALSE)[, c(1L, 2L)]
      }
    }, error = function(e) NULL)
    if (!is.null(mm_raw) && ncol(mm_raw) >= 2) {
      colnames(mm_raw) <- c("cluster_id", "protein_id")
      mm_cluster_df <- unique(mm_raw[mm_raw$cluster_id %in% cluster_ids & !is.na(mm_raw$protein_id) & mm_raw$protein_id != "", , drop = FALSE])
    }
  }

  # --- 5. Token-Level Waterfall Resolution (Method 1 -> Method 2 -> Method 3) ---
  empty_tier <- data.frame(
    cluster_id = character(0), target = character(0),
    match_tier = character(0), tier_rank = integer(0),
    stringsAsFactors = FALSE
  )

  # Method 1: exact token match
  m1_join <- unique(merge(cluster_tok_df[, c("cluster_id", "token", "pair_key")], m_tok_df[, c("target", "token")], by = "token"))
  m1_keys <- unique(m1_join$pair_key)
  m1_hits <- if (nrow(m1_join) > 0) {
    unique(data.frame(
      cluster_id = m1_join$cluster_id,
      target = m1_join$target,
      match_tier = "Method1_Direct_LocusSymbol",
      tier_rank = 1L,
      stringsAsFactors = FALSE
    ))
  } else {
    empty_tier
  }

  # Method 2: evaluate only unmapped (cluster_id, token) pairs via AureoWiki
  rem1_tok <- cluster_tok_df[!cluster_tok_df$pair_key %in% m1_keys, , drop = FALSE]
  m2_keys <- character(0)
  m2_hits <- empty_tier
  if (nrow(rem1_tok) > 0 && nrow(pan_to_m_df) > 0) {
    c_to_pan <- unique(merge(rem1_tok[, c("cluster_id", "token", "pair_key")], pan_tok_df[, c("pan_id", "token")], by = "token")[, c("cluster_id", "pan_id", "pair_key")])
    if (nrow(c_to_pan) > 0) {
      m2_join <- unique(merge(c_to_pan, pan_to_m_df, by = "pan_id"))
      if (nrow(m2_join) > 0) {
        m2_keys <- unique(m2_join$pair_key)
        m2_hits <- unique(data.frame(
          cluster_id = m2_join$cluster_id,
          target = m2_join$target,
          match_tier = "Method2_AureoWiki_Ortholog",
          tier_rank = 2L,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Method 3: evaluate remaining (cluster_id, token) pairs + MMseqs2 alignments
  rem2_tok <- rem1_tok[!rem1_tok$pair_key %in% m2_keys, c("cluster_id", "token_low"), drop = FALSE]
  m3_hits <- empty_tier
  t3_parts <- list()
  if (nrow(rem2_tok) > 0) {
    t3_parts[[1]] <- unique(merge(rem2_tok, m_tok_for_pan[, c("target", "token_low")], by = "token_low")[, c("cluster_id", "target")])
    if (nrow(pan_to_m_df) > 0) {
      c_to_pan_ci <- unique(merge(rem2_tok, pan_tok_df[, c("pan_id", "token_low")], by = "token_low")[, c("cluster_id", "pan_id")])
      if (nrow(c_to_pan_ci) > 0) {
        t3_parts[[length(t3_parts) + 1]] <- unique(merge(c_to_pan_ci, pan_to_m_df, by = "pan_id")[, c("cluster_id", "target")])
      }
    }
  }
  if (nrow(mm_cluster_df) > 0) {
    mm_keys <- paste0(mm_cluster_df$cluster_id, "\t", mm_cluster_df$protein_id)
    mm_sub <- mm_cluster_df[!mm_keys %in% c(m1_keys, m2_keys), , drop = FALSE]
    if (nrow(mm_sub) > 0) {
      t3_parts[[length(t3_parts) + 1]] <- unique(merge(mm_sub, m_tok_for_pan[, c("target", "token")], by.x = "protein_id", by.y = "token")[, c("cluster_id", "target")])
      if (nrow(pan_to_m_df) > 0) {
        mm_pan <- unique(merge(mm_sub, pan_tok_df[, c("pan_id", "token")], by.x = "protein_id", by.y = "token")[, c("cluster_id", "pan_id")])
        if (nrow(mm_pan) > 0) {
          t3_parts[[length(t3_parts) + 1]] <- unique(merge(mm_pan, pan_to_m_df, by = "pan_id")[, c("cluster_id", "target")])
        }
      }
    }
  }
  t3_df <- if (length(t3_parts) > 0) unique(do.call(rbind, t3_parts)) else data.frame(cluster_id = character(0), target = character(0), stringsAsFactors = FALSE)
  if (nrow(t3_df) > 0) {
    m3_hits <- data.frame(
      cluster_id = t3_df$cluster_id,
      target = t3_df$target,
      match_tier = "Method3_MMseqs2_GAF_Fallback",
      tier_rank = 3L,
      stringsAsFactors = FALSE
    )
  }

  # Keep highest-priority tier per unique (cluster_id, target) pair
  raw_map <- rbind(m1_hits, m2_hits, m3_hits)
  if (nrow(raw_map) > 0) {
    raw_map <- raw_map[order(raw_map$cluster_id, raw_map$target, raw_map$tier_rank), , drop = FALSE]
    raw_map <- raw_map[!duplicated(raw_map[, c("cluster_id", "target")]), , drop = FALSE]
  }
  mapping_df <- raw_map[, c("cluster_id", "target", "match_tier"), drop = FALSE]
  rownames(mapping_df) <- NULL

  if (nrow(raw_map) > 0) {
    agg_tier <- aggregate(match_tier ~ cluster_id, data = raw_map, FUN = function(x) {
      ord <- c("Method1_Direct_LocusSymbol", "Method2_AureoWiki_Ortholog", "Method3_MMseqs2_GAF_Fallback")
      u <- unique(x)
      paste(c(intersect(ord, u), setdiff(u, ord)), collapse = " + ")
    })
    agg_tgt <- aggregate(target ~ cluster_id, data = raw_map, FUN = function(x) paste(unique(x), collapse = "/"))
    cnt_tgt <- aggregate(target ~ cluster_id, data = raw_map, FUN = function(x) length(unique(x)))
    colnames(agg_tgt)[2] <- "mapped_targets"
    colnames(cnt_tgt)[2] <- "n_targets"
    mapped_sum <- Reduce(function(a, b) merge(a, b, by = "cluster_id"), list(agg_tier, agg_tgt, cnt_tgt))
  } else {
    mapped_sum <- data.frame(cluster_id = character(0), match_tier = character(0), mapped_targets = character(0), n_targets = integer(0), stringsAsFactors = FALSE)
  }

  cluster_summary <- merge(data.frame(cluster_id = cluster_ids, stringsAsFactors = FALSE), mapped_sum, by = "cluster_id", all.x = TRUE, sort = FALSE)
  cluster_summary$n_targets[is.na(cluster_summary$n_targets)] <- 0L
  cluster_summary$match_tier[is.na(cluster_summary$match_tier)] <- "Unmapped"

  n_total <- length(cluster_ids)
  n_mapped <- sum(cluster_summary$n_targets > 0)
  n_multi <- sum(cluster_summary$n_targets > 1)
  pct_mapped <- if (n_total > 0) 100 * n_mapped / n_total else 0
  pct_multi <- if (n_total > 0) 100 * n_multi / n_total else 0
  n_m_covered <- length(unique(mapping_df$target))
  pct_m_covered <- if (length(m_gene_ids) > 0) 100 * n_m_covered / length(m_gene_ids) else 0
  tier_counts <- table(cluster_summary$match_tier)

  if (verbose) {
    cat(sprintf(
      "  ✓ Cluster-to-M mapping: %d / %d cluster_ids (%.2f%%) mapped to %d / %d M-matrix genes (%.2f%%).\n",
      n_mapped, n_total, pct_mapped, n_m_covered, length(m_gene_ids), pct_m_covered
    ))
    cat(sprintf(
      "    Multi-target clusters (>1 USA300HOU_RS*): %d / %d (%.2f%%)\n",
      n_multi, n_total, pct_multi
    ))
    cat("    Tier breakdown -> ",
        paste(sprintf("%s: %d", names(tier_counts), as.integer(tier_counts)), collapse = " | "),
        "\n", sep = "")
  }

  list(
    mapping_df = mapping_df,
    cluster_summary = cluster_summary,
    stats = list(
      n_total_clusters = n_total,
      n_mapped_clusters = n_mapped,
      pct_mapped_clusters = pct_mapped,
      n_multi_target_clusters = n_multi,
      pct_multi_target_clusters = pct_multi,
      n_total_m_genes = length(m_gene_ids),
      n_mapped_m_genes = n_m_covered,
      pct_mapped_m_genes = pct_m_covered,
      tier_counts = as.list(tier_counts)
    )
  )
}
