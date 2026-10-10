import requests
import json
import time
import os

with open('../unique_protein_ids.txt') as f:
    input_ids = f.read().split('\t')

print(f"Total input IDs: {len(input_ids)}")

chunk_size = 200
url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
results = {}

for i in range(0, len(input_ids), chunk_size):
    chunk = input_ids[i:i+chunk_size]
    params = {'db': 'protein', 'id': ','.join(chunk), 'retmode': 'json'}
    try:
        res = requests.post(url, data=params).json()
        if 'result' in res and 'uids' in res['result']:
            for uid in res['result']['uids']:
                if uid == 'error':
                    continue
                record = res['result'][uid]
                acc_ver = record.get('accessionversion', '')
                if not acc_ver:
                    continue
                
                active = acc_ver
                if record.get('replacedby'):
                    active = record['replacedby']
                elif record.get('identicalto') and isinstance(record['identicalto'], dict) and 'acc' in record['identicalto']:
                     active = record['identicalto']['acc']
                
                results[acc_ver] = active
        
    except Exception as e:
        print(f"Error in batch {i}: {e}")
    time.sleep(0.34) # ~3 requests per second limit

# Now find redundancies
active_to_originals = {}
for orig, active in results.items():
    active_to_originals.setdefault(active, []).append(orig)

redundant_groups = {active: origs for active, origs in active_to_originals.items() if len(origs) > 1}

with open('redundant_mappings.txt', 'w') as f:
    for active, origs in redundant_groups.items():
        f.write(f"Active Protein {active} represents duplicates: {', '.join(origs)}\n")

print(f"Found {len(redundant_groups)} groups of redundant IDs.")
redundant_count = sum(len(origs) - 1 for origs in redundant_groups.values())
print(f"Total redundant IDs mapped: {redundant_count}")
