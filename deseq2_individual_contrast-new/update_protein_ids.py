import os
import pandas as pd

def add_new_protein_ids(result_dir):
    mapping_file = os.path.join(result_dir, 'protein_id_mapping.csv')
    if not os.path.exists(mapping_file):
        print(f"Error: {mapping_file} not found. Please run create_id_mapping.py first.")
        return

    # Load mapping dictionary
    print(f"Loading ID mapping from {mapping_file}...")
    try:
        mapping_df = pd.read_csv(mapping_file)
        # Convert to dictionary for fast lookup
        id_mapping = dict(zip(mapping_df['original_protein_id'], mapping_df['fasta_protein_id']))
        print(f"Loaded {len(id_mapping)} mappings.")
    except Exception as e:
        print(f"Error loading mapping file: {e}")
        return

    target_files = []
    for subfolder in ['raw', 'microarray']:
        sub_dir = os.path.join(result_dir, subfolder)
        if os.path.exists(sub_dir):
            for filename in os.listdir(sub_dir):
                if filename.endswith(('.csv', '.xlsx', '.xls')):
                    target_files.append(os.path.join(subfolder, filename))
    
    if not target_files:
        print("No data files found in raw/ or microarray/ subdirectories.")
        return

    for rel_path in target_files:
        if os.path.basename(rel_path).startswith(('add_cluster_info_summary', 'add_gene_info_summary')):
            continue
            
        full_path = os.path.join(result_dir, rel_path)
        dataset_name = rel_path.replace(os.sep, '/')
        
        print(f"Processing {dataset_name}...")
        try:
            if full_path.endswith('.csv'):
                df = pd.read_csv(full_path)
            else:
                df = pd.read_excel(full_path)
            
            if 'protein_id' not in df.columns:
                print(f"  [WARNING] 'protein_id' column not found. Skipping file.")
                continue
                
            # Count valid protein IDs originally
            has_protein_id = df['protein_id'].notna() & (df['protein_id'].astype(str).str.strip() != "")
            total_protein_ids = has_protein_id.sum()
            
            # Map the IDs
            def resolve_id(old_id):
                if pd.isna(old_id) or str(old_id).strip() == "":
                    return ""
                
                # Fetch mapped value
                new_val = id_mapping.get(str(old_id).strip(), "")
                
                # Blank out if it was an error/dropped record
                if new_val in ["Not found in FASTA", "Invalid/Not Found"]:
                    return ""
                return new_val

            df['new_protein_id'] = df['protein_id'].apply(resolve_id)
            
            # Count how many successfully got a new ID
            has_new_id = df['new_protein_id'] != ""
            added_new_ids = has_new_id.sum()
            
            # Save it back (overwrite)
            if full_path.endswith('.csv'):
                df.to_csv(full_path, index=False)
            else:
                df.to_excel(full_path, index=False)
            
            print(f"  [SUCCESS] Added 'new_protein_id'. Mapped {added_new_ids} / {total_protein_ids} valid protein IDs.")
            
        except Exception as e:
            print(f"  [ERROR] Failed to process {dataset_name}: {e}")

if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    add_new_protein_ids(current_dir)
