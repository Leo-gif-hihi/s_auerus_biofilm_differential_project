import os
import glob
from collections import defaultdict
import urllib.parse

def check_overlaps(annote_dir):
    # Mapping of tag -> set of filenames
    tag_to_files = defaultdict(set)
    # Mapping of filename -> set of tags (for reporting)
    file_to_tags = defaultdict(set)

    gff_files = glob.glob(os.path.join(annote_dir, "*.gff3"))
    
    print(f"Checking {len(gff_files)} .gff3 files for overlapping locus tags...")
    
    for gff_path in gff_files:
        filename = os.path.basename(gff_path)
        with open(gff_path, 'r', encoding='utf-8') as f:
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
                
                # Extract both new and old locus tags
                tags = set()
                if 'locus_tag' in attributes:
                    tags.update([t.strip() for t in attributes['locus_tag'].split(',') if t.strip()])
                if 'old_locus_tag' in attributes:
                    tags.update([t.strip() for t in attributes['old_locus_tag'].split(',') if t.strip()])
                    
                for tag in tags:
                    tag_to_files[tag].add(filename)
                    file_to_tags[filename].add(tag)
                    
    # Find overlapping tags
    overlapping_tags = {tag: files for tag, files in tag_to_files.items() if len(files) > 1}
    
    if not overlapping_tags:
        print("\nNo overlapping locus tags found between different .gff3 files.")
        return
        
    print(f"\nFound {len(overlapping_tags)} locus tags that appear in multiple annotation files.")
    
    # Group by the specific file combinations to make it readable
    combination_to_tags = defaultdict(list)
    for tag, files in overlapping_tags.items():
        combo = tuple(sorted(list(files)))
        combination_to_tags[combo].append(tag)
        
    for combo, tags in sorted(combination_to_tags.items(), key=lambda x: len(x[1]), reverse=True):
        print(f"\n==================================================")
        print(f"Overlap between: {', '.join(combo)}")
        print(f"Total shared tags: {len(tags)}")
        
        # Print up to 20 examples to not flood the console
        sample_tags = tags[:20]
        print(f"Examples: {', '.join(sample_tags)}")
        if len(tags) > 20:
            print(f"... and {len(tags) - 20} more")

if __name__ == "__main__":
    annote_directory = os.path.join(os.path.dirname(__file__), 'annote')
    check_overlaps(annote_directory)
