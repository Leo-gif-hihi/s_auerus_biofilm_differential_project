import requests
import json
import time
import pandas as pd

# 1. Load original IDs
with open('unique_protein_ids.txt', 'r', encoding='utf-8') as f:
    original_ids = f.read().split('\t')

# 2. Load FASTA IDs
with open('all_protein_sequences.fasta', 'r', encoding='utf-8') as f:
    fasta_ids = set([line.split()[0][1:] for line in f if line.startswith('>')])

# 3. Direct matches
mapping = {}
missing_ids = []

for oid in original_ids:
    if oid in fasta_ids:
        mapping[oid] = oid
    else:
        missing_ids.append(oid)

print(f"Direct matches: {len(mapping)}")
print(f"Missing IDs to resolve: {len(missing_ids)}")

# 4. Resolve missing via NCBI
url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
chunk_size = 200

# We need to submit batches of missing IDs
for i in range(0, len(missing_ids), chunk_size):
    chunk = missing_ids[i:i+chunk_size]
    params = {'db': 'protein', 'id': ','.join(chunk), 'retmode': 'json'}
    
    try:
        res = requests.post(url, data=params).json()
        if 'result' in res and 'uids' in res['result']:
            for uid in res['result']['uids']:
                if uid == 'error': continue
                record = res['result'][uid]
                
                acc_ver = record.get('accessionversion', '')
                if not acc_ver: continue
                
                # Check for replacements
                active = acc_ver
                if record.get('replacedby'):
                    active = record['replacedby']
                elif record.get('identicalto') and isinstance(record['identicalto'], dict) and 'acc' in record['identicalto']:
                     active = record['identicalto']['acc']
                
                # We need to map the INPUT ID to the ACTIVE ID
                # Since we don't strictly know which input mapped to this record (if accessionversion doesn't match),
                # but accessionversion is almost always the input ID:
                if acc_ver in missing_ids:
                    mapping[acc_ver] = active if active in fasta_ids else 'Not found in FASTA'
                
                # Alternative: Some old records might omit version
                acc_no_ver = acc_ver.split('.')[0]
                if acc_no_ver in missing_ids and acc_ver not in missing_ids:
                    mapping[acc_no_ver] = active if active in fasta_ids else 'Not found in FASTA'

    except Exception as e:
        print(f"Error in batch {i}: {e}")
    
    time.sleep(0.34)

# 5. Output mapping
# For any missing ID that didn't get resolved (e.g., completely invalid), mark them
for oid in missing_ids:
    if oid not in mapping:
        mapping[oid] = 'Invalid/Not Found'

# Save to CSV
df = pd.DataFrame(list(mapping.items()), columns=['original_protein_id', 'fasta_protein_id'])
df.to_csv('protein_id_mapping.csv', index=False)

print(f"Successfully saved mapping for {len(df)} IDs to protein_id_mapping.csv")
