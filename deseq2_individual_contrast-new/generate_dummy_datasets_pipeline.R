set.seed(42)

base_dir <- "simulated_pipeline_data"
raw_dir <- file.path(base_dir, "raw")
sample_dir <- file.path(base_dir, "sample_info-all")

if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)
if (!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)

num_datasets <- 8
num_genes <- 2500
prop_sig <- 0.5

# Generate generic cluster IDs and locus tags
cluster_ids <- paste0("cluster_", 1:num_genes)
locus_tags <- paste0("SAOUHSC_", sprintf("%05d", 1:num_genes))

# We need an annot.csv that has "GSE ID" and "Modified in biofilm formation compare to control (description)"
annot_df <- data.frame(
  GSE_ID = character(),
  Direction = character(),
  stringsAsFactors = FALSE
)

for (i in 1:num_datasets) {
  # The pipeline script uses `str_extract(..., "GSE[0-9]+")` to find dataset IDs
  gse_id <- paste0("GSE99999", i)
  
  # Randomly select 50% of genes to be significant in THIS specific dataset
  sig_genes_idx <- sample(1:num_genes, size = floor(num_genes * prop_sig))
  
  sample_names <- paste0(rep(c("Control", "Treatment"), each=3), "_rep", 1:3)
  
  counts <- matrix(0, nrow = num_genes, ncol = 6)
  colnames(counts) <- sample_names
  
  for (j in 1:num_genes) {
    base_lambda <- sample(50:1000, 1) 
    
    if (j %in% sig_genes_idx) {
      is_up <- sample(c(TRUE, FALSE), 1)
      effect_size <- sample(3:6, 1) 
      
      if (is_up) {
        counts[j, 1:3] <- rpois(3, lambda = base_lambda)
        counts[j, 4:6] <- rpois(3, lambda = base_lambda * effect_size)
      } else {
        counts[j, 1:3] <- rpois(3, lambda = base_lambda * effect_size)
        counts[j, 4:6] <- rpois(3, lambda = base_lambda)
      }
    } else {
      counts[j, ] <- rpois(6, lambda = base_lambda)
    }
  }
  
  # 1. Format raw counts exactly as the pipeline expects
  raw_df <- data.frame(
    cluster_id = cluster_ids,
    gene = locus_tags,
    counts,
    check.names = FALSE
  )
  write.csv(raw_df, file.path(raw_dir, paste0(gse_id, "_counts.csv")), row.names = FALSE)
  
  # 2. Format sample info exactly as the pipeline expects
  si_df <- data.frame(
    Sample = sample_names,
    Group = rep(c("Control", "Treatment"), each=3),
    check.names = FALSE
  )
  write.csv(si_df, file.path(sample_dir, paste0(gse_id, "_sample_info.csv")), row.names = FALSE)
  
  # 3. Add to annotation dataframe
  direction <- sample(c("increased", "decreased"), 1)
  annot_df <- rbind(annot_df, data.frame(
    GSE_ID = gse_id,
    Direction = direction,
    stringsAsFactors = FALSE
  ))
  
  # Save the list of true DE genes for verification
  true_de_genes <- data.frame(
    cluster_id = cluster_ids[sig_genes_idx], 
    gene = locus_tags[sig_genes_idx]
  )
  write.csv(true_de_genes, file.path(base_dir, paste0(gse_id, "_true_DE_genes.csv")), row.names = FALSE)
}

# Write out the annot.csv with correct column names
colnames(annot_df) <- c("GSE ID", "Modified in biofilm formation compare to control (description)")
write.csv(annot_df, file.path(base_dir, "annot.csv"), row.names = FALSE, quote = FALSE)

cat("Successfully generated pipeline-compatible dummy data in:", base_dir, "\n")
