run_ora_meta_analysis <- function(gaf_file, lfc_list, pval_list, cts, cluster_annotation, comparison_gse_map, min_ds, out_dir, ora_padj = 0.05, ora_lfc = 1.0) {
  cat("\n=== RUNNING OVER-REPRESENTATION META-ANALYSIS ===\n")
  if (!is.null(gaf_file) && file.exists(gaf_file)) {
    cat("Loading GAF file:", gaf_file, "\n")
    
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)
    }
    if (!requireNamespace("GO.db", quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install("GO.db", update = FALSE, ask = FALSE)
    }
    
    suppressPackageStartupMessages(library(clusterProfiler))
    suppressPackageStartupMessages(library(GO.db))
    
    gaf <- read.delim(gaf_file, comment.char = "!", header = FALSE, stringsAsFactors = FALSE, quote="")
    term2gene <- unique(data.frame(term = gaf$V5, gene = gaf$V3, stringsAsFactors = FALSE))
    
    annotated_genes <- unique(gaf$V3)
    
    cluster_id_to_symbol <- sapply(cluster_annotation$All_Genes_Locus_Tags, function(tags) {
      if (is.na(tags) || tags == "") return(NA_character_)
      tag_list <- strsplit(as.character(tags), "/")[[1]]
      matches <- tag_list[tag_list %in% annotated_genes]
      if (length(matches) > 0) {
        return(matches[1])
      } else {
        return(NA_character_)
      }
    })
    names(cluster_id_to_symbol) <- cluster_annotation$cluster_id
    
    goterms <- as.data.frame(GO.db::GOTERM)
    term2name <- unique(data.frame(term = goterms$go_id, name = goterms$Term, stringsAsFactors = FALSE))
    
    ora_logor_list <- list()
    ora_se_list <- list()
    ora_comparisons <- c()
    ora_gse_map <- c()
    
    for (comp_name in names(lfc_list)) {
      gse_id <- comparison_gse_map[comp_name]
      
      pvals <- pval_list[[comp_name]]
      lfcs <- lfc_list[[comp_name]]
      
      valid_mask <- !is.na(pvals) & !is.na(lfcs)
      universe_clusters <- rownames(cts)[valid_mask]
      
      if (length(universe_clusters) == 0) next
      
      dataset_padj <- p.adjust(pvals[valid_mask], method = "BH")
      de_mask <- dataset_padj < ora_padj & abs(lfcs[valid_mask]) >= ora_lfc
      de_clusters <- universe_clusters[de_mask]
      
      if (length(de_clusters) == 0) {
        cat("  - No DE genes for", comp_name, "- skipping ORA.\n")
        next
      }
      
      universe_symbols <- unique(cluster_id_to_symbol[universe_clusters])
      universe_symbols <- universe_symbols[!is.na(universe_symbols) & universe_symbols != ""]
      
      de_symbols <- unique(cluster_id_to_symbol[de_clusters])
      de_symbols <- de_symbols[!is.na(de_symbols) & de_symbols != ""]
      
      if (length(de_symbols) == 0) next
      
      ora_res <- tryCatch({
        enricher(
          gene = de_symbols,
          universe = universe_symbols,
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          minGSSize = 5,
          maxGSSize = 500,
          TERM2GENE = term2gene,
          TERM2NAME = term2name
        )
      }, error = function(e) NULL)
      
      if (!is.null(ora_res) && nrow(ora_res@result) > 0) {
        res_df <- ora_res@result
        
        parse_num <- function(x) as.numeric(sapply(strsplit(as.character(x), "/"), `[`, 1))
        parse_den <- function(x) as.numeric(sapply(strsplit(as.character(x), "/"), `[`, 2))
        
        a <- parse_num(res_df$GeneRatio)
        de_size <- parse_den(res_df$GeneRatio)
        b <- de_size - a
        
        bg_num <- parse_num(res_df$BgRatio)
        bg_den <- parse_den(res_df$BgRatio)
        c_val <- bg_num - a
        d <- bg_den - de_size - c_val
        
        a_adj <- a + 0.5
        b_adj <- b + 0.5
        c_adj <- c_val + 0.5
        d_adj <- d + 0.5
        
        logOR <- log((a_adj * d_adj) / (b_adj * c_adj))
        se_logOR <- sqrt(1/a_adj + 1/b_adj + 1/c_adj + 1/d_adj)
        
        ora_logor_list[[comp_name]] <- setNames(logOR, res_df$ID)
        ora_se_list[[comp_name]] <- setNames(se_logOR, res_df$ID)
        ora_comparisons <- c(ora_comparisons, comp_name)
        ora_gse_map[comp_name] <- gse_id
      }
    }
    
    if (length(ora_logor_list) > 0) {
      all_go_ids <- unique(unlist(lapply(ora_logor_list, names)))
      
      ora_yi <- matrix(NA, nrow = length(all_go_ids), ncol = length(ora_logor_list), dimnames = list(all_go_ids, ora_comparisons))
      ora_vi <- matrix(NA, nrow = length(all_go_ids), ncol = length(ora_se_list), dimnames = list(all_go_ids, ora_comparisons))
      
      for (comp_name in ora_comparisons) {
        ora_yi[names(ora_logor_list[[comp_name]]), comp_name] <- ora_logor_list[[comp_name]]
        ora_vi[names(ora_se_list[[comp_name]]), comp_name] <- ora_se_list[[comp_name]]^2
      }
      
      cat("  Meta-analyzing", length(all_go_ids), "GO terms across", length(ora_comparisons), "comparisons...\n")
      
      ora_meta_res <- list()
      for (i in seq_along(all_go_ids)) {
        go_id <- all_go_ids[i]
        yi_go <- ora_yi[i, ]
        vi_go <- ora_vi[i, ]
        
        valid_go <- !is.na(yi_go) & !is.na(vi_go)
        valid_gse_go <- ora_gse_map[valid_go]
        valid_comp_go <- ora_comparisons[valid_go]
        
        if (sum(valid_go) >= max(2, min_ds)) {
          res_re <- tryCatch(
            metafor::rma.mv(yi = yi_go[valid_go], V = vi_go[valid_go], random = ~ 1 | valid_gse_go / valid_comp_go, method = "REML"),
            error = function(e) {
              tryCatch(
                metafor::rma.mv(yi = yi_go[valid_go], V = vi_go[valid_go], random = ~ 1 | valid_gse_go, method = "REML"),
                error = function(e2) NULL
              )
            }
          )
          
          res_fe <- tryCatch(
            metafor::rma(yi = yi_go[valid_go], vi = vi_go[valid_go], method = "FE"),
            error = function(e) NULL
          )
          
          if (!is.null(res_re) && !is.null(res_fe)) {
            ora_meta_res[[go_id]] <- data.frame(
              ID = go_id,
              Description = term2name$name[match(go_id, term2name$term)],
              RE_logOR = res_re$beta[1],
              RE_SE = res_re$se[1],
              RE_pval = res_re$pval[1],
              FE_logOR = res_fe$beta[1],
              FE_SE = res_fe$se[1],
              FE_pval = res_fe$pval[1],
              N_comparisons = sum(valid_go),
              N_datasets = length(unique(valid_gse_go)),
              stringsAsFactors = FALSE
            )
          }
        }
      }
      
      if (length(ora_meta_res) > 0) {
        ora_meta_df <- do.call(rbind, ora_meta_res)
        ora_meta_df <- ora_meta_df[order(ora_meta_df$RE_pval), ]
        write.csv(ora_meta_df, file.path(out_dir, "ORA_Meta_Analysis_Results.csv"), row.names = FALSE)
        cat("  - Saved ORA Meta-Analysis Results to", file.path(out_dir, "ORA_Meta_Analysis_Results.csv"), "\n")
      } else {
        cat("  - Not enough valid datasets/comparisons to run meta-analysis on any GO term.\n")
      }
    } else {
      cat("  - No ORA results could be computed for any comparison.\n")
    }
  } else {
    cat("  - GAF file not found or not provided. Skipping ORA meta-analysis.\n")
  }
}
