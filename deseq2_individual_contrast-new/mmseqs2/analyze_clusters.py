import csv

def get_sequence_lengths(fasta_file):
    lengths = {}
    with open(fasta_file, 'r') as f:
        curr_id = None
        curr_len = 0
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if curr_id:
                    lengths[curr_id] = curr_len
                curr_id = line[1:].split()[0]
                curr_len = 0
            else:
                curr_len += len(line)
        if curr_id:
            lengths[curr_id] = curr_len
    return lengths

def analyze_clusters(input_file, fasta_file, output_file):
    unique_clusters = set()
    cluster_sizes = {}
    target_to_queries = {}
    
    min_identity = float('inf')
    min_identity_pairs = []
    
    min_qcov = float('inf')
    min_qcov_pairs = []

    # Bins for identity
    id_bins = {
        '1.000': 0,
        '0.990 - 0.999': 0,
        '0.950 - 0.989': 0,
        '0.900 - 0.949': 0,
        '0.800 - 0.899': 0,
        '0.700 - 0.799': 0,
        '< 0.700': 0
    }

    diff_len_count = 0
    same_len_count = 0
    diff_len_examples = []
    
    seq_lengths = get_sequence_lengths(fasta_file)

    with open(input_file, 'r', newline='') as f:
        fieldnames = ["query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "qcov", "tcov"]
        reader = csv.DictReader(f, delimiter='\t', fieldnames=fieldnames)
        for row in reader:
            query = row['query']
            target = row['target']
            
            # The cluster is typically represented by the query sequence
            unique_clusters.add(query)
            cluster_sizes[query] = cluster_sizes.get(query, 0) + 1
            
            if target not in target_to_queries:
                target_to_queries[target] = set()
            target_to_queries[target].add(query)
            
            fident = float(row['fident'])
            qcov = float(row['qcov'])
            
            if query != target:
                # Identity binning
                if fident >= 1.0:
                    id_bins['1.000'] += 1
                elif fident >= 0.99:
                    id_bins['0.990 - 0.999'] += 1
                elif fident >= 0.95:
                    id_bins['0.950 - 0.989'] += 1
                elif fident >= 0.90:
                    id_bins['0.900 - 0.949'] += 1
                elif fident >= 0.80:
                    id_bins['0.800 - 0.899'] += 1
                elif fident >= 0.70:
                    id_bins['0.700 - 0.799'] += 1
                else:
                    id_bins['< 0.700'] += 1
                    
                # Length calculation
                qlen = seq_lengths.get(query, 0)
                tlen = seq_lengths.get(target, 0)
                
                if qlen != tlen:
                    diff_len_count += 1
                    if len(diff_len_examples) < 5:
                        diff_len_examples.append({'query': query, 'qlen': qlen, 'target': target, 'tlen': tlen, 'fident': fident})
                else:
                    same_len_count += 1

            # 4. Lowest identity
            if fident < min_identity:
                min_identity = fident
                min_identity_pairs = [{'query': query, 'target': target, 'fident': fident}]
            elif fident == min_identity:
                min_identity_pairs.append({'query': query, 'target': target, 'fident': fident})
                
            # 4. Lowest qcov
            if qcov < min_qcov:
                min_qcov = qcov
                min_qcov_pairs = [{'query': query, 'target': target, 'qcov': qcov}]
            elif qcov == min_qcov:
                min_qcov_pairs.append({'query': query, 'target': target, 'qcov': qcov})

    # Write results to output file
    with open(output_file, 'w') as out:
        out.write("--- Cluster Alignment Analysis ---\n\n")
        
        # 1. Unique clusters
        out.write(f"1. How many unique clusters are there?\n")
        out.write(f"Answer: {len(unique_clusters)}\n\n")
        
        # 2. Sequence Identity Breakdown
        out.write("2. Sequence Identity Breakdown:\n")
        out.write("Sequence Identity Range\tAlignment Pairs\n")
        out.write(f"1.000 (100% identical)\t{id_bins['1.000']}\n")
        out.write(f"0.990 - 0.999\t{id_bins['0.990 - 0.999']}\n")
        out.write(f"0.950 - 0.989\t{id_bins['0.950 - 0.989']}\n")
        out.write(f"0.900 - 0.949\t{id_bins['0.900 - 0.949']}\n")
        out.write(f"0.800 - 0.899\t{id_bins['0.800 - 0.899']}\n")
        out.write(f"0.700 - 0.799\t{id_bins['0.700 - 0.799']}\n")
        out.write(f"< 0.700\t{id_bins['< 0.700']}\n\n")
        
        # 3. Different sequence lengths
        total_clustered_pairs = diff_len_count + same_len_count
        same_len_perc = (same_len_count / total_clustered_pairs * 100) if total_clustered_pairs > 0 else 0
        out.write("3. Length differences:\n")
        out.write(f"There are proteins of different lengths grouped into the same cluster in your data (specifically, {diff_len_count} pairs have different lengths), but ~{same_len_perc:.0f}% of your clustered pairs have the exact same length.\n\n")
        if diff_len_examples:
            out.write("Here are real examples from your cluster_alignments.tsv where lengths differ:\n\n")
            first_q = diff_len_examples[0]['query']
            first_qlen = diff_len_examples[0]['qlen']
            out.write(f"Representative {first_q} (length {first_qlen} aa) clusters with:\n")
            for ex in diff_len_examples:
                if ex['query'] == first_q:
                    out.write(f"{ex['target']} (length {ex['tlen']} aa)\n")
        out.write("\n")
        
        # 4. Lowest identity
        out.write(f"4a. What is the lowest value of identity?\n")
        out.write(f"Answer: {min_identity}\n")
        out.write("Which comparison hold the lowest value?\n")
        for item in min_identity_pairs:
            out.write(f"Query: {item['query']}, Target: {item['target']}, fident: {item['fident']}\n")
        out.write("\n")
        
        # 4. Lowest qcov
        out.write(f"4b. What is the lowest value of qcov?\n")
        out.write(f"Answer: {min_qcov}\n")
        out.write("Which comparison hold the lowest value?\n")
        for item in min_qcov_pairs:
            out.write(f"Query: {item['query']}, Target: {item['target']}, qcov: {item['qcov']}\n")
        out.write("\n")
        
        # 5. Single protein clusters
        single_protein_clusters = sum(1 for size in cluster_sizes.values() if size == 1)
        out.write(f"5. How many clusters containing only one protein?\n")
        out.write(f"Answer: {single_protein_clusters}\n\n")
        
        # 6. Targets mapped to multiple clusters
        multi_cluster_targets = sum(1 for queries in target_to_queries.values() if len(queries) > 1)
        out.write(f"6. How many target proteins map to multiple clusters?\n")
        out.write(f"Answer: {multi_cluster_targets}\n\n")
        
    print(f"Analysis complete! Results saved to {output_file}")

if __name__ == '__main__':
    input_tsv = 'cluster_alignments.tsv'
    fasta_file = 'all_protein_sequences.fasta'
    output_txt = 'cluster_analysis_results.txt'
    analyze_clusters(input_tsv, fasta_file, output_txt)
