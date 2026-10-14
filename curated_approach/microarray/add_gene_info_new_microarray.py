import os
import argparse
import pandas as pd
import re
import glob
import urllib.parse
import requests
import time
import sys

class TeeLogger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a", encoding='utf-8')

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()

log_path = os.path.join(os.path.dirname(__file__), "add_gene_info.log")
sys.stdout = TeeLogger(log_path)

def is_valid_val(val):
    if val is None or pd.isna(val):
        return False
    s = str(val).strip()
    return s not in ['', '-', '—', '–', '―', '−', '--', '---', '.', 'nan', 'NaN', 'none', 'None', 'null', 'NULL', 'NA', 'N/A', '#N/A']

def fetch_aureowiki_info(gene_symbol):
    """
    Queries Aureowiki API for a gene symbol and extracts product, protein accession ID,
    and equivalent identifiers. Returns a dictionary with the extracted data.
    """
    api_url = "https://aureowiki.med.uni-greifswald.de/api.php"
    
    # Search for the page title
    search_params = {
        "action": "query",
        "list": "search",
        "srsearch": gene_symbol,
        "utf8": "1",
        "format": "json",
        "srlimit": 20
    }
    
    try:
        response = requests.get(api_url, params=search_params, headers={"User-Agent": "AureoWikiQueryBot/1.0"})
        data = response.json()
        search_results = data.get("query", {}).get("search", [])
        
        if not search_results:
            return None
            
        for i in range(min(5, len(search_results))):
            page_title = search_results[i]["title"]
            original_page_title = page_title
            saupan_identifiers = []
            
            if "SAUPAN" in page_title:
                # Fetch SAUPAN page to extract ortholog link
                parse_params = {
                    "action": "parse",
                    "page": page_title,
                    "redirects": "1",
                    "prop": "text",
                    "format": "json"
                }
                response = requests.get(api_url, params=parse_params, headers={"User-Agent": "AureoWikiQueryBot/1.0"})
                data = response.json()
                if "parse" in data and "text" in data["parse"]:
                    saupan_html = data["parse"]["text"]["*"]
                    
                    import html
                    matches = re.findall(r'<div class="[^"]*ab_orthotable_value_cell[^"]*"[^>]*>(.*?)</div>', saupan_html, re.DOTALL)
                    for m in matches:
                        raw_text = re.sub(r'<[^>]+>', ' ', m)
                        raw_text = html.unescape(raw_text).strip()
                        if raw_text and raw_text != '—':
                            tokens = re.split(r'[\s\(\)]+', raw_text)
                            for t in tokens:
                                if t:
                                    saupan_identifiers.append(t)
                                    
                    links = re.findall(r'ab_orthotable_value_cell.*?href="/([^"]+)"', saupan_html, re.DOTALL)
                    
                    valid_link = None
                    for link in links:
                        if "index.php" not in link and ":" not in link:
                            valid_link = link
                            break
                            
                    if valid_link:
                        page_title = valid_link
                    else:
                        continue
                else:
                    continue
        
        # Fetch rendered HTML
            # Fetch rendered HTML
            parse_params = {
                "action": "parse",
                "page": page_title,
                "redirects": "1",
                "prop": "text",
                "format": "json"
            }
            
            response = requests.get(api_url, params=parse_params, headers={"User-Agent": "AureoWikiQueryBot/1.0"})
            data = response.json()
            
            if "parse" not in data or "text" not in data["parse"]:
                continue
                
            html_content = data["parse"]["text"]["*"]
            
            current_query = gene_symbol
            discontinued_match = re.search(r'(?i)discontinued.*?<a[^>]+href="/([^"]+)"', html_content, re.DOTALL)
            if discontinued_match:
                new_gene = discontinued_match.group(1)
                print(f"Gene {current_query} was discontinued on page {page_title}. Replacing with {new_gene}...")
                current_query = new_gene
                page_title = new_gene
                
                parse_params["page"] = page_title
                response = requests.get(api_url, params=parse_params, headers={"User-Agent": "AureoWikiQueryBot/1.0"})
                data = response.json()
                if "parse" not in data or "text" not in data["parse"]:
                    continue
                html_content = data["parse"]["text"]["*"]
            
            result = {
                'locus tag': None,
                'symbol': None,
                'pan gene symbol': None,
                'synonym': None,
                'alternate name': None,
                'product': None,
                'protein RefSeq': None,
                'type': None,
                'organism': None,
                'pan locus tag': None,
                'protein sequence': None,
                'essential': None
            }
            
            import html
            
            lis = re.findall(r'<li[^>]*>(.*?)</li>', html_content, re.DOTALL | re.IGNORECASE)
            for li in lis:
                text = re.sub(r'<[^>]+>', ' ', li)
                text = re.sub(r'\s+', ' ', text).strip()
                text = html.unescape(text)
                text_lower = text.lower()
                
                if text_lower.startswith("locus tag:"):
                    val = text.split(":", 1)[-1].strip()
                    if is_valid_val(val):
                        # might contain "[new locus tag: ...]" so we take the first word
                        result['locus tag'] = val.split()[0]
                elif text_lower.startswith("symbol:"):
                    val = text.split(":", 1)[-1].strip()
                    if is_valid_val(val): result['symbol'] = val
                elif text_lower.startswith("pan gene symbol"):
                    val = text.split(":", 1)[-1].strip()
                    if is_valid_val(val): result['pan gene symbol'] = val
                elif text_lower.startswith("synonym:"):
                    val = text.split(":", 1)[-1].strip()
                    if is_valid_val(val): result['synonym'] = val
                elif text_lower.startswith("alternate name:"):
                    val = text.split(":", 1)[-1].strip()
                    if is_valid_val(val): result['alternate name'] = val
                elif text_lower.startswith("product:"):
                    val = text.split(":", 1)[-1].strip()
                    if is_valid_val(val): result['product'] = val
                elif text_lower.startswith("type:"):
                    val = text.split(":", 1)[-1].strip()
                    if is_valid_val(val): result['type'] = val
                elif text_lower.startswith("refseq:"):
                    val = text.split(":", 1)[-1].strip()
                    if is_valid_val(val): result['protein RefSeq'] = val.split()[0]
                elif text_lower.startswith("organism:"):
                    val = text.split(":", 1)[-1].strip()
                    if is_valid_val(val): result['organism'] = val
                elif text_lower.startswith("pan locus tag"):
                    val = text.split(":", 1)[-1].strip()
                    if is_valid_val(val): result['pan locus tag'] = val
                elif text_lower.startswith("essential:"):
                    val = text.split(":", 1)[-1].strip()
                    if is_valid_val(val): result['essential'] = val
                        
            seq_match = re.search(r'(?i)id="Protein_sequence".*?<(?:ul|ol)>.*?<li[^>]*>(.*?)</li>', html_content, re.DOTALL)
            if seq_match:
                seq = html.unescape(re.sub(r'<[^>]+>', '', seq_match.group(1))).strip()
                if seq: result['protein sequence'] = seq
                
            if not any(result.values()):
                continue
                
            def strip_final_digit(s):
                if not s: return ""
                return re.sub(r'\d$', '', str(s)).strip().lower()
                
            def validate_aureowiki_match(query, result_dict, p_title, saupan_ids):
                if not result_dict: return False
                q_lower = str(query).strip().lower()
                q_stripped = strip_final_digit(q_lower)
                
                targets = [p_title]
                targets.extend(saupan_ids)
                for k in ['locus tag', 'pan gene symbol', 'symbol', 'synonym', 'alternate name']:
                    if result_dict.get(k):
                        val_str = str(result_dict[k])
                        for v in val_str.replace(',', '/').split('/'):
                            targets.append(v.strip())
                        
                for t in targets:
                    if not t: continue
                    t_str = str(t).strip().lower()
                    t_stripped = strip_final_digit(t_str)
                    if q_lower == t_str or q_stripped == t_stripped:
                        return True
                return False
                
            if not validate_aureowiki_match(current_query, result, original_page_title, saupan_identifiers):
                print(f"Validation failed for {current_query} on result {i+1} against AureoWiki results. Trying next...")
                continue
                
            return result
            
        print(f"All validation attempts failed for {gene_symbol}. Treating as Not Found.")
        return None
        
    except Exception as e:
        print(f"Error fetching from Aureowiki for {gene_symbol}: {e}")
        return None


def merge_gene_info(existing, new_info):
    if existing['priority'] < new_info['priority']:
        return {
            'symbol': existing['symbol'] if existing['symbol'] else new_info['symbol'],
            'product': existing['product'] if existing['product'] else new_info['product'],
            'protein_id': existing['protein_id'] if existing['protein_id'] else new_info['protein_id'],
            'is_pseudogene': existing['is_pseudogene'],
            'priority': existing['priority'],
            'protein_sequence': existing.get('protein_sequence') if existing.get('protein_sequence') else new_info.get('protein_sequence')
        }
    elif existing['priority'] > new_info['priority']:
        return {
            'symbol': new_info['symbol'] if new_info['symbol'] else existing['symbol'],
            'product': new_info['product'] if new_info['product'] else existing['product'],
            'protein_id': new_info['protein_id'] if new_info['protein_id'] else existing['protein_id'],
            'is_pseudogene': new_info['is_pseudogene'],
            'priority': new_info['priority'],
            'protein_sequence': new_info.get('protein_sequence') if new_info.get('protein_sequence') else existing.get('protein_sequence')
        }
    else:
        return {
            'symbol': existing['symbol'] if existing['symbol'] else new_info['symbol'],
            'product': existing['product'] if existing['product'] else new_info['product'],
            'protein_id': existing['protein_id'] if existing['protein_id'] else new_info['protein_id'],
            'is_pseudogene': existing['is_pseudogene'] or new_info['is_pseudogene'],
            'priority': existing['priority'],
            'protein_sequence': existing.get('protein_sequence') if existing.get('protein_sequence') else new_info.get('protein_sequence')
        }

def parse_gff(gff_path, mapping=None, priority=3):
    """
    Parses a GFF file and updates the dictionary mapping identifiers (locus_tag, old_locus_tag, gene symbol)
    to a dictionary containing {symbol, product, protein_id, is_pseudogene, priority}.
    """
    if mapping is None:
        mapping = {}
    
    print(f"Parsing {os.path.basename(gff_path)} (Priority {priority})...")
    
    locus_agg = {}
    
    try:
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                feature_type = parts[2]
                if feature_type not in ['CDS', 'gene', 'pseudogene', 'mRNA', 'tRNA', 'rRNA']:
                    continue
                
                attributes_str = parts[8]
                attributes = {}
                for attr in attributes_str.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key] = urllib.parse.unquote(value)
                
                locus_tag = attributes.get('locus_tag')
                if not locus_tag:
                    locus_tag = attributes.get('ID')
                    
                if not locus_tag:
                    continue
                    
                if locus_tag not in locus_agg:
                    locus_agg[locus_tag] = {
                        'symbol': None,
                        'product': None,
                        'protein_id': None,
                        'is_pseudogene': False,
                        'aliases': set()
                    }
                
                agg = locus_agg[locus_tag]
                
                for key in ['old_locus_tag', 'gene', 'locus_tag']:
                    if attributes.get(key):
                        for val in attributes[key].split(','):
                            alias_val = val.strip()
                            if alias_val:
                                # Fix NCBI annotation error where N315 locus tags (SAxxxx) 
                                # are incorrectly listed as old locus tags in COL.gff3
                                if 'COL.gff3' in gff_path and re.match(r'^SA\d{4}$', alias_val):
                                    continue
                                agg['aliases'].add(alias_val)
                    
                if attributes.get('gene') and not agg['symbol'] and is_valid_val(attributes.get('gene')):
                    agg['symbol'] = attributes.get('gene').strip()
                if attributes.get('product') and is_valid_val(attributes.get('product')):
                    agg['product'] = attributes.get('product').strip()
                if attributes.get('protein_id') and is_valid_val(attributes.get('protein_id')):
                    agg['protein_id'] = attributes.get('protein_id').strip()
                    
                is_pseudogene = (
                    feature_type == 'pseudogene' or
                    attributes.get('pseudo') == 'true' or
                    attributes.get('gene_biotype') == 'pseudogene'
                )
                if is_pseudogene:
                    agg['is_pseudogene'] = True
                    
        for locus_tag, agg in locus_agg.items():
            info = {
                'symbol': agg['symbol'],
                'product': agg['product'],
                'protein_id': agg['protein_id'],
                'is_pseudogene': agg['is_pseudogene'],
                'priority': priority,
                'protein_sequence': None
            }
            
            if locus_tag not in mapping:
                mapping[locus_tag] = info.copy()
            else:
                mapping[locus_tag] = merge_gene_info(mapping[locus_tag], info)
                
            for alias in agg['aliases']:
                if alias not in mapping:
                    mapping[alias] = info.copy()
                else:
                    mapping[alias] = merge_gene_info(mapping[alias], info)
                
    except Exception as e:
        print(f"Error parsing {gff_path}: {e}")
        
    return mapping

def parse_excel(excel_path, mapping=None, priority=4, queried_set=None):
    if mapping is None:
        mapping = {}
        
    print(f"Parsing {os.path.basename(excel_path)} (Priority {priority})...")
    
    try:
        df = pd.read_excel(excel_path)
        
        for _, row in df.iterrows():
            locus_tag = str(row.get('locus tag', '')).strip()
            old_locus_tag = str(row.get('old locus tag', '')).strip()
            pan_symbol = str(row.get('pan gene symbol', '')).strip()
            symbol = str(row.get('symbol', '')).strip()
            synonym = str(row.get('synonym', '')).strip()
            queried_tag = str(row.get('queried_tag', '')).strip()
            
            product = str(row.get('product', '')).strip()
            if not is_valid_val(product): product = None
            
            protein_id = str(row.get('protein RefSeq', '')).strip()
            if not is_valid_val(protein_id): protein_id = None

            protein_sequence = str(row.get('protein sequence', '')).strip()
            if not is_valid_val(protein_sequence): protein_sequence = None
            
            primary_symbol = None
            if is_valid_val(pan_symbol):
                primary_symbol = pan_symbol
            elif is_valid_val(symbol):
                primary_symbol = symbol
                
            is_pseudogene = False
            if 'type' in df.columns:
                type_val = str(row.get('type', '')).strip().lower()
                if 'pseudogene' in type_val:
                    is_pseudogene = True
                    
            aliases = set()
            for val in [locus_tag, old_locus_tag, pan_symbol, symbol, synonym, queried_tag]:
                if is_valid_val(val):
                    aliases.add(val)
                    
            if is_valid_val(queried_tag) and queried_set is not None:
                queried_set.add(queried_tag)
                
            if not any([primary_symbol, product, protein_id, is_pseudogene]):
                continue
                
            if not aliases:
                continue
                
            info = {
                'symbol': primary_symbol,
                'product': product,
                'protein_id': protein_id,
                'is_pseudogene': is_pseudogene,
                'priority': priority,
                'protein_sequence': protein_sequence
            }
            
            for alias in aliases:
                if alias not in mapping:
                    mapping[alias] = info.copy()
                else:
                    mapping[alias] = merge_gene_info(mapping[alias], info)
                    
    except Exception as e:
        print(f"Error parsing {excel_path}: {e}")
        
    return mapping

def process_datasets(skip_aureowiki_tags=None, target_file=None):
    if skip_aureowiki_tags is None:
        skip_aureowiki_tags = []
        
    # Ensure we use the script's location as the root directory to find 'annote' reliably
    base_dir = os.path.dirname(os.path.abspath(__file__))
    result_dir = base_dir
    annote_dir = os.path.join(base_dir, 'annote')
    
    if not os.path.exists(result_dir) or not os.path.exists(annote_dir):
        print("Result or annote directory not found.")
        return
    
    summary_stats = []
    
    # Store fasta sequences across all datasets
    fasta_entries = {}
    
    normalized_log = []
    seen_normalized = set()
    
    def clean_tags(gene_id_str, seen_norm, norm_log):
        if pd.isna(gene_id_str):
            return []
        tags_list = [t.strip() for t in re.split(r'[,;/]+', str(gene_id_str)) if t.strip()]
        tags_list = [t[5:] if t.startswith('gene-') else t for t in tags_list]
        tags_list = [re.sub(r'\s*\(.*?\)', '', t).strip() for t in tags_list]
        
        new_tags = []
        for t in tags_list:
            match = re.match(r'^([a-zA-Z]+)_(\d{1,2})$', t)
            if match:
                norm_tag = f"{match.group(1)}{match.group(2)}"
                new_tags.append(norm_tag)
                log_key = f"{t}->{norm_tag}"
                if log_key not in seen_norm:
                    seen_norm.add(log_key)
                    norm_log.append({
                        'Original Tag': t,
                        'Normalized Tag': norm_tag,
                        'Found in Mapping': 'Pending search'
                    })
            else:
                new_tags.append(t)
                
        explicit_map = {
            'kdpC(SCCmec)' : 'kdpC'
        }
        explicit_map_lower = {k.lower(): v for k, v in explicit_map.items()}
        final_tags = []
        for t in new_tags:
            t_lower = t.lower()
            if t_lower in explicit_map_lower:
                final_tags.append(explicit_map_lower[t_lower])
            else:
                final_tags.append(t)
                
        return final_tags
    
    aureowiki_new_results = []
    queried_aureowiki = set()
    
    print("Pre-parsing all annotation files into cache...")
    annote_cache = {}
    for filename in os.listdir(annote_dir):
        file_path = os.path.join(annote_dir, filename)
        if filename.endswith('.gff') or filename.endswith('.gff3'):
            mapping = parse_gff(file_path, priority=99)
            annote_cache[filename] = {
                'tags': set(mapping.keys()),
                'mapping': mapping,
                'type': 'gff'
            }
        elif filename.endswith('.xlsx') and filename.startswith('GeneSpecific'):
            mapping = parse_excel(file_path, priority=99, queried_set=queried_aureowiki)
            annote_cache[filename] = {
                'tags': set(mapping.keys()),
                'mapping': mapping,
                'type': 'excel'
            }

    target_files = []
    if target_file:
        target_files.append(os.path.abspath(target_file))
    else:
        platforms = ['Affymatrix', 'Agilent', 'NimboGen']
        for platform in platforms:
            platform_dir = os.path.join(base_dir, platform)
            if os.path.isdir(platform_dir):
                for root, dirs, files in os.walk(platform_dir):
                    for file in files:
                        # Look for Differential Expression CSV files
                        if "_DE_" in file and file.endswith('.csv'):
                            target_files.append(os.path.join(root, file))

    for file_path in target_files:
        filename = os.path.basename(file_path)
        try:
            rel_path = os.path.relpath(file_path, base_dir)
        except ValueError:
            rel_path = file_path
            
        print(f"Processing {rel_path}...")
        
        try:
            if filename.endswith('.csv'):
                df = pd.read_csv(file_path, encoding='utf-8-sig')
            else:
                df = pd.read_excel(file_path)
                
            # Clean up artifacts from previous re-runs
            if 'is_pseudogene' in df.columns:
                df.drop(columns=['gene', 'product', 'protein_id', 'is_pseudogene'], inplace=True, errors='ignore')
                
            corrupted_cols = [c for c in df.columns if re.match(r'^(locus_tag|gene)\.\d+$', str(c))]
            if corrupted_cols:
                df.drop(columns=corrupted_cols, inplace=True, errors='ignore')
            
            possible_gene_cols = ['locus_tag', 'gene_id', 'gene', 'gene_name', 'final_gene_name']
            gene_col = next((col for col in possible_gene_cols if col in df.columns), None)
            
            if not gene_col:
                gene_col = df.columns[0]
                print(f"  Warning: Expected gene column not found in {filename}. Using first column '{gene_col}' as gene identifiers.")
            
            dataset_locus_tags = set()
            for gene_id in df[gene_col]:
                tags = clean_tags(gene_id, seen_normalized, normalized_log)
                for t in tags:
                    if len(t) >= 6 and re.match(r'^[a-zA-Z0-9_\-]+$', t):
                        dataset_locus_tags.add(t)
                        
            overlap_counts = {}
            for annote_file, cache_data in annote_cache.items():
                if annote_file == 'GeneSpecific_Aureowiki.xlsx':
                    continue
                overlap = dataset_locus_tags.intersection(cache_data['tags'])
                overlap_counts[annote_file] = len(overlap)
                
            max_overlap = max(overlap_counts.values()) if overlap_counts else 0
            
            prioritized_files = []
            if max_overlap > 0:
                threshold = max(max_overlap * 0.2, min(max_overlap, 5))
                for annote_file, count in overlap_counts.items():
                    if count >= threshold:
                        prioritized_files.append(annote_file)
                    
            print(f"  Prioritized annotation files for {filename}: {prioritized_files if prioritized_files else 'None'}")
            
            global_mapping = {}
            def merge_mapping(dest_mapping, src_mapping, priority):
                for k, info in src_mapping.items():
                    new_info = info.copy()
                    new_info['priority'] = priority
                    if k not in dest_mapping:
                        dest_mapping[k] = new_info
                    else:
                        existing = dest_mapping[k]
                        dest_mapping[k] = merge_gene_info(existing, new_info)
            
            processed_files = set()
            for annote_file in prioritized_files:
                if annote_file in annote_cache:
                    cache_data = annote_cache[annote_file]
                    prio = 1 if cache_data['type'] == 'gff' else 2
                    merge_mapping(global_mapping, cache_data['mapping'], prio)
                    processed_files.add(annote_file)
            
            for annote_file in sorted(annote_cache.keys()):
                if annote_file not in processed_files:
                    cache_data = annote_cache[annote_file]
                    prio = 3 if cache_data['type'] == 'gff' else 4
                    merge_mapping(global_mapping, cache_data['mapping'], prio)
            
            new_genes = []
            products = []
            protein_ids = []
            is_pseudogenes = []
            
            stats = {
                'Dataset': rel_path,
                'Total Genes': len(df),
                'Product Added': 0,
                'Protein ID Added': 0,
                'Pseudogene Added': 0
            }
            
            for gene_id in df[gene_col]:
                gene_id_str = str(gene_id).strip()
                tags = clean_tags(gene_id, seen_normalized, normalized_log)
                
                matched_infos = [
                    (t, global_mapping[t]) 
                    for t in tags if t in global_mapping
                ]
                
                if not matched_infos:
                    for t in tags:
                        should_skip = False
                        for skip_str in skip_aureowiki_tags:
                            if skip_str in t:
                                should_skip = True
                                break
                                
                        if should_skip:
                            continue
                            
                        if t not in queried_aureowiki:
                            print(f"    Querying AureoWiki for {t}...")
                            time.sleep(1.5)
                            aureo_data = fetch_aureowiki_info(t)
                            queried_aureowiki.add(t)
                            
                            record = aureo_data.copy() if aureo_data else {}
                            record['queried_tag'] = t
                            aureowiki_new_results.append(record)
                            
                            # Save immediately to prevent data loss
                            cache_file = os.path.join(annote_dir, 'GeneSpecific_Aureowiki.xlsx')
                            try:
                                df_new = pd.DataFrame([record])
                                if os.path.exists(cache_file):
                                    df_exist = pd.read_excel(cache_file)
                                    df_combined = pd.concat([df_exist, df_new]).drop_duplicates()
                                else:
                                    df_combined = df_new
                                df_combined.to_excel(cache_file, index=False)
                            except Exception as e:
                                print(f"    Warning: Could not save AureoWiki cache for {t}: {e}")
                            
                            if aureo_data:
                                if aureo_data.get('protein RefSeq') or aureo_data.get('protein sequence'):
                                    type_val = str(aureo_data.get('type')).strip().lower()
                                    is_pseudo = ('pseudogene' in type_val)
                                    info = {
                                        'symbol': aureo_data.get('symbol') or aureo_data.get('pan gene symbol'),
                                        'product': aureo_data.get('product'),
                                        'protein_id': aureo_data.get('protein RefSeq'),
                                        'is_pseudogene': is_pseudo,
                                        'priority': 5,
                                        'protein_sequence': aureo_data.get('protein sequence')
                                    }
                                    matched_infos.append((t, info))
                                    global_mapping[t] = info
                                    break
                                else:
                                    eq_ids = []
                                    for key in ['locus tag', 'symbol', 'pan gene symbol', 'synonym', 'alternate name']:
                                        if aureo_data.get(key):
                                            val_str = str(aureo_data[key])
                                            for v in val_str.replace(',', '/').split('/'):
                                                eq_ids.append(v.strip())
                                            
                                    found_eq = False
                                    for eq_id in eq_ids:
                                        if eq_id in global_mapping:
                                            matched_infos.append((eq_id, global_mapping[eq_id]))
                                            global_mapping[t] = global_mapping[eq_id]
                                            found_eq = True
                                            break
                                    if found_eq:
                                        break
                                        

                if not matched_infos:
                    for t in tags:
                        # 1. Rule: Remove numbers (arcC1 -> arcC)
                        match = re.match(r'^([a-zA-Z]+)_?\d{1,2}$', t)
                        if match:
                            norm_tag = match.group(1)
                            found = norm_tag in global_mapping
                            
                            if t not in seen_normalized:
                                seen_normalized.add(t)
                                normalized_log.append({
                                    'Original Tag': t,
                                    'Normalized Tag': norm_tag,
                                    'Found in Mapping': found
                                })
                                
                            if found:
                                matched_infos.append((norm_tag, global_mapping[norm_tag]))
                        
                        # 2. Rule: Add numbers (arcC -> arcC1, arcC_1, etc.)
                        if t.isalpha():
                            for i in range(1, 11):
                                for suffix in [str(i), f"_{i}"]:
                                    norm_tag = f"{t}{suffix}"
                                    found = norm_tag in global_mapping
                                    
                                    if found:
                                        log_key = f"{t}->{norm_tag}"
                                        if log_key not in seen_normalized:
                                            seen_normalized.add(log_key)
                                            normalized_log.append({
                                                'Original Tag': t,
                                                'Normalized Tag': norm_tag,
                                                'Found in Mapping': found
                                            })
                                        matched_infos.append((norm_tag, global_mapping[norm_tag]))
                                
                matched_infos.sort(key=lambda x: x[1]['priority'])
                
                if matched_infos:
                    primary_tag, primary_info = matched_infos[0]
                    
                    symbol = next((info['symbol'] for _, info in matched_infos if is_valid_val(info.get('symbol'))), None)
                    product = next((info['product'] for _, info in matched_infos if is_valid_val(info.get('product'))), None)
                    protein_id = next((info['protein_id'] for _, info in matched_infos if is_valid_val(info.get('protein_id'))), None)
                    protein_sequence = next((info.get('protein_sequence') for _, info in matched_infos if is_valid_val(info.get('protein_sequence'))), None)
                    best_priority = primary_info['priority']
                    is_pseudogene = any(info.get('is_pseudogene', False) for _, info in matched_infos if info['priority'] == best_priority)
                    
                    new_gene = symbol if is_valid_val(symbol) else primary_tag
                    if not is_valid_val(new_gene):
                        new_gene = tags[0] if (tags and is_valid_val(tags[0])) else gene_id_str
                    new_genes.append(new_gene)
                    
                    products.append(product)
                    if product:
                        stats['Product Added'] += 1
                        
                    if not protein_id and protein_sequence:
                        unique_name = f"{primary_tag}_protein"
                        protein_id = unique_name
                        fasta_entries[unique_name] = protein_sequence

                    protein_ids.append(protein_id)
                    if protein_id:
                        stats['Protein ID Added'] += 1
                        
                    is_pseudogenes.append(is_pseudogene)
                    if is_pseudogene:
                        stats['Pseudogene Added'] += 1
                else:
                    fallback_tag = tags[0] if (tags and is_valid_val(tags[0])) else gene_id_str
                    new_genes.append(fallback_tag)
                    products.append(None)
                    protein_ids.append(None)
                    is_pseudogenes.append(None)
            
            if gene_col == 'gene':
                if 'locus_tag' in df.columns:
                    df.rename(columns={'gene': 'original_gene'}, inplace=True)
                    gene_col = 'original_gene'
                else:
                    df.rename(columns={'gene': 'locus_tag'}, inplace=True)
                    gene_col = 'locus_tag'
                
            col_idx = df.columns.tolist().index(gene_col) + 1
            df.insert(col_idx, 'gene', new_genes)
            
            df['product'] = products
            df['protein_id'] = protein_ids
            df['is_pseudogene'] = is_pseudogenes
            
            total = stats['Total Genes']
            if total > 0:
                stats['Product Added (%)'] = round((stats['Product Added'] / total) * 100, 2)
                stats['Protein ID Added (%)'] = round((stats['Protein ID Added'] / total) * 100, 2)
                stats['Pseudogene Added (%)'] = round((stats['Pseudogene Added'] / total) * 100, 2)
            else:
                stats['Product Added (%)'] = 0.0
                stats['Protein ID Added (%)'] = 0.0
                stats['Pseudogene Added (%)'] = 0.0

            if filename.endswith('.csv'):
                df.to_csv(file_path, index=False, encoding='utf-8-sig')
            else:
                df.to_excel(file_path, index=False)
                
            summary_stats.append(stats)
            print(f"  Updated {rel_path}")
            
        except Exception as e:
            print(f"  Error processing {rel_path}: {e}")

    if summary_stats:
        summary_df = pd.DataFrame(summary_stats)
        summary_path = os.path.join(result_dir, 'add_gene_info_summary.xlsx')
        try:
            if os.path.exists(summary_path):
                existing_df = pd.read_excel(summary_path)
                summary_df = pd.concat([existing_df, summary_df], ignore_index=True)
                summary_df.drop_duplicates(subset=['Dataset'], keep='last', inplace=True)
            summary_df.to_excel(summary_path, index=False)
            print(f"\nSummary saved to {summary_path}")
        except Exception as e:
            print(f"\nCould not save to Excel ({e}). Saving to CSV instead...")
            summary_path_csv = os.path.join(result_dir, 'add_gene_info_summary.csv')
            if os.path.exists(summary_path_csv):
                existing_df = pd.read_csv(summary_path_csv)
                summary_df = pd.concat([existing_df, summary_df], ignore_index=True)
                summary_df.drop_duplicates(subset=['Dataset'], keep='last', inplace=True)
            summary_df.to_csv(summary_path_csv, index=False)
            print(f"Summary saved to {summary_path_csv}")

    if normalized_log:
        norm_df = pd.DataFrame(normalized_log)
        norm_path = os.path.join(result_dir, 'normalized_tags_log.csv')
        norm_df.to_csv(norm_path, index=False)
        print(f"Normalized tags log saved to {norm_path}")
        
    if fasta_entries:
        fasta_path = os.path.join(base_dir, 'missing_protein_id.fasta')
        try:
            with open(fasta_path, 'w') as f:
                for unique_name, seq in fasta_entries.items():
                    f.write(f">{unique_name}\n")
                    # Optionally wrap sequence to 80 chars, but a single line is fine too
                    f.write(f"{seq}\n")
            print(f"Fasta sequences exported to {fasta_path}")
        except Exception as e:
            print(f"Error saving FASTA file: {e}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add gene info to microarray datasets')
    parser.add_argument('--skip-aureowiki-tags', nargs='*', default=["UAMS-1_", "SA564_", "SKAG_", "DA471_", "SAUSA300s", "Predicted_", "Novel", "sRNA", "EBG0000", "AVU41_"], 
                        help='Strings to skip when querying AureoWiki (e.g., UAMS-1_)')
    parser.add_argument('--file', type=str, help='Specific file to process', default=None)
    args = parser.parse_args()
    
    process_datasets(skip_aureowiki_tags=args.skip_aureowiki_tags, target_file=args.file)
