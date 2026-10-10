import os
import glob
import pandas as pd

def add_cluster_information(result_dir, tsv_file):
    print(f"Loading cluster mapping from {tsv_file}...")
    try:
        # We only need the first two columns: 'query' (cluster_id) and 'target' (protein_id)
        # MMseqs2 default TSV output has no headers
        cluster_df = pd.read_csv(tsv_file, sep='\t', header=None, usecols=[0, 1], names=['query', 'target'])
    except Exception as e:
        print(f"Error loading TSV file: {e}")
        return

    # Create a mapping dictionary from target (protein_id) to query (cluster_id)
    # Dropping duplicates to ensure a 1:1 mapping (target -> query)
    cluster_df = cluster_df.drop_duplicates(subset=['target'])
    protein_to_cluster = dict(zip(cluster_df['target'], cluster_df['query']))
    print(f"Loaded mapping for {len(protein_to_cluster)} unique proteins.")

    target_files = []
    for subfolder in ['raw', 'microarray']:
        sub_dir = os.path.join(result_dir, subfolder)
        if os.path.exists(sub_dir):
            for filename in os.listdir(sub_dir):
                if filename.endswith('.csv') or filename.endswith('.xlsx') or filename.endswith('.xls'):
                    target_files.append(os.path.join(subfolder, filename))
    
    if not target_files:
        print("No data files found in the directory.")
        return

    summary_stats = []

    for rel_path in target_files:
        if os.path.basename(rel_path).startswith('add_cluster_info_summary'):
            continue
            
        full_path = os.path.join(result_dir, rel_path)
        dataset_name = rel_path.replace(os.sep, '/')
        
        print(f"Processing {dataset_name}...")
        try:
            if full_path.endswith('.csv'):
                df = pd.read_csv(full_path)
            else:
                df = pd.read_excel(full_path)
            
            # Check if the required column exists
            if 'new_protein_id' not in df.columns:
                print(f"  [WARNING] 'new_protein_id' column not found. Skipping file.")
                continue
                
            # Count valid protein IDs (not NaN and not empty string)
            has_protein_id = df['new_protein_id'].notna() & (df['new_protein_id'].astype(str).str.strip() != "")
            total_protein_ids = has_protein_id.sum()
            
            # Map the cluster_id using the protein_to_cluster dictionary
            # Any new_protein_id not found in the dictionary will be left empty
            df['cluster_id'] = df['new_protein_id'].map(protein_to_cluster).fillna("")
            
            # Count successfully mapped cluster IDs
            has_cluster_id = df['cluster_id'] != ""
            added_cluster_ids = has_cluster_id.sum()
            
            # Calculate percentage
            percentage = (added_cluster_ids / total_protein_ids * 100) if total_protein_ids > 0 else 0.0
            
            # Overwrite the original file
            if full_path.endswith('.csv'):
                df.to_csv(full_path, index=False)
            else:
                df.to_excel(full_path, index=False)
            
            print(f"  [STAT] Cluster IDs added: {added_cluster_ids} / {total_protein_ids} protein IDs ({percentage:.2f}%)")
            print(f"  [SUCCESS] Updated {dataset_name} with 'cluster_id' column.")
            
            # Collect stats for summary
            summary_stats.append({
                'Dataset': dataset_name,
                'Total Protein IDs': total_protein_ids,
                'Cluster IDs Added': added_cluster_ids,
                'Percentage (%)': round(percentage, 2)
            })
            
        except Exception as e:
            print(f"  [ERROR] Failed to process {dataset_name}: {e}")

    # Save summary report if there are any processed files
    if summary_stats:
        summary_df = pd.DataFrame(summary_stats)
        summary_csv_path = os.path.join(result_dir, "add_cluster_info_summary.csv")
        summary_df.to_csv(summary_csv_path, index=False)
        print(f"\nSaved summary report to {summary_csv_path}")

if __name__ == "__main__":
    # The script is designed to be run from inside the 'result' directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Try the exact name first, and fallback to the one without 's' just in case
    tsv_path = os.path.join(current_dir, "mmseqs2", "cluster_alignments.tsv")
    if not os.path.exists(tsv_path):
        tsv_path = os.path.join(current_dir, "mmseq2", "cluster_alignments.tsv")
        
    add_cluster_information(current_dir, tsv_path)
