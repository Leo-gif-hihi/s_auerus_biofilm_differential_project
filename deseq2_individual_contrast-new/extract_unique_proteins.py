import os
import pandas as pd

def extract_unique_proteins(directory):
    unique_proteins = set()
    data_files = []
    for subfolder in ['raw', 'microarray']:
        sub_dir = os.path.join(directory, subfolder)
        if os.path.exists(sub_dir):
            for f in os.listdir(sub_dir):
                if f.endswith('.csv') or f.endswith('.xlsx') or f.endswith('.xls'):
                    data_files.append(os.path.join(subfolder, f))
    
    for file in data_files:
        filepath = os.path.join(directory, file)
        try:
            # Read only the header first to see if protein_id is in it
            if filepath.endswith('.csv'):
                df_header = pd.read_csv(filepath, nrows=0)
            else:
                df_header = pd.read_excel(filepath, nrows=0)
                
            if 'protein_id' in df_header.columns:
                # Read the actual column
                if filepath.endswith('.csv'):
                    df = pd.read_csv(filepath, usecols=['protein_id'])
                else:
                    df = pd.read_excel(filepath, usecols=['protein_id'])
                    
                # Drop NA and empty strings
                valid_proteins = df['protein_id'].dropna().astype(str).str.strip()
                valid_proteins = valid_proteins[(valid_proteins != '') & (valid_proteins != 'nan') & (valid_proteins != 'NA') & (valid_proteins != 'None')]
                unique_proteins.update(valid_proteins)
        except Exception as e:
            print(f"Error processing {file}: {e}")
            
    # Write to file separated by tab
    output_path = os.path.join(directory, 'unique_protein_ids.txt')
    with open(output_path, 'w') as f:
        f.write('\t'.join(sorted(unique_proteins)))
        
    print(f"Number of unique protein_id: {len(unique_proteins)}")

if __name__ == '__main__':
    extract_unique_proteins('.')
