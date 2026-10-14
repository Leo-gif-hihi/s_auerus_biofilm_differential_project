# Microarray Gene Annotation Pipeline (`add_gene_info_new_microarray.py`)

This README details the logic, data flow, and specific edge-case handling implemented in `add_gene_info_new_microarray.py`. This script is designed to standardize gene symbols, fill in missing gene products, and assign RefSeq protein IDs (`WP_...`) to microarray differential expression datasets across different *Staphylococcus aureus* strains.

---

## 1. Core Logic Overview

For each microarray dataset (CSV/Excel files found in the `Affymatrix`, `Agilent`, and `NimboGen` directories), the script performs the following steps:
1. **Reads the dataset** and identifies the column containing gene identifiers (e.g., `locus_tag`, `gene_id`, `gene`).
2. **Cleans and normalizes** the gene identifiers.
3. **Prioritizes annotation databases** by calculating the overlap between the dataset's locus tags and all available `.gff3` and `.xlsx` annotation files in the `annote/` directory.
4. **Builds a `global_mapping`** dictionary merging data from all annotation files, respecting the calculated priorities.
5. **Maps each gene** in the dataset to the best available `symbol`, `product`, and `protein_id` (RefSeq).
6. **Falls back to AureoWiki API** for any unresolved tags, caching the results locally to avoid redundant API calls.
7. **Overwrites the dataset** with the newly appended `gene`, `product`, `protein_id`, and `is_pseudogene` columns.

---

## 2. Tag Normalization (`clean_tags`)

Microarray datasets often have messy gene names. The `clean_tags` function standardizes them:
- Splits tags delimited by commas, semicolons, or slashes (e.g., `SA0184; SACOL0176`).
- Strips the `gene-` prefix if it exists.
- Removes parenthetical text (except for explicitly hardcoded exceptions).
- Normalizes underscores before digits (e.g., `arcC_1` is transformed to `arcC1`).
- Contains a hardcoded explicit map (e.g., `kdpC(SCCmec)` -> `kdpC`).

---

## 3. Annotation File Parsing & Prioritization

The script relies on reference files located in the `annote/` directory:
- **`.gff3` files**: Official NCBI annotations.
- **`.xlsx` files**: Curated gene specific information tables (`GeneSpecificInformation_*.xlsx`).

### Prioritization Algorithm
Because microarray datasets span different strains (e.g., USA300, COL, N315), the script automatically identifies the closest matching annotation files:
- It checks the overlap between the normalized tags in the dataset and the tags present in each annotation file.
- Files exceeding an overlap threshold (20% of the maximum overlap) are marked as **Prioritized**.

When building the `global_mapping` dictionary, it resolves conflicts based on the following priority hierarchy (Lower number = Higher priority):
- **Priority 1**: Prioritized `.gff3` files
- **Priority 2**: Prioritized `.xlsx` files
- **Priority 3**: Non-prioritized `.gff3` files
- **Priority 4**: Non-prioritized `.xlsx` files
- **Priority 5**: AureoWiki API fallback cache

### Data Merging (`merge_gene_info`)
When the script encounters the same locus tag or alias in multiple files, it merges them. If an attribute (like `protein_id`) is missing in a Priority 1 file, it will inherit it from a Priority 2 or 3 file, ensuring no data is lost.

### Important Edge Case Fix (`parse_gff`)
There is a known NCBI annotation error where the `COL.gff3` file incorrectly lists N315 locus tags (which follow the format `SA\d{4}`, e.g., `SA0184`) as `old_locus_tag` aliases for completely different genes (e.g., `gisA`).
- **The Fix**: The script explicitly filters out any alias matching the `^SA\d{4}$` regex exclusively when parsing `COL.gff3`. This prevents cross-contamination where N315 genes (like `mupG`) were having their symbols merged with COL protein IDs (like `WP_000067349.1`).

---

## 4. Alias Resolution and Fallbacks

If a normalized tag is found in the `global_mapping` dictionary, the script assigns the data. If it is **not found**, the script attempts multiple fallbacks:

### A. AureoWiki API Query (`fetch_aureowiki_info`)
- The script queries the `aureowiki.med.uni-greifswald.de` API.
- It parses the HTML payload of the wiki page to extract `locus tag`, `symbol`, `product`, `protein RefSeq`, and `type`.
- It tracks orthologous pages (e.g., redirecting to SAUPAN pan-genome pages or following "discontinued" gene notices).
- It strictly validates the search result against the queried tag to ensure a false positive wasn't returned.
- Successful queries are cached in `annote/GeneSpecific_Aureowiki.xlsx` so they are never queried twice.

### B. Fuzzy Fallbacks (Digit Stripping/Adding)
If AureoWiki also fails, the script assumes the gene is an ortholog variant:
1. **Rule 1 (Remove numbers)**: If the tag is `arcC1`, it tries looking up `arcC` in the `global_mapping`.
2. **Rule 2 (Add numbers)**: If the tag is `arcC`, it iterates from 1 to 10 and checks if `arcC1`, `arcC_1`, `arcC2`, etc., exist in the `global_mapping`.

### C. Global Mapping Search Keys
The search against `global_mapping` involves matching dataset queries against specific keys stored in the mapping.
- **Dataset Query Terms Searched**:
  - Direct normalized tags from the dataset's gene identifier column.
  - AureoWiki equivalent identifiers (Fallback A).
  - Fuzzy ortholog variants like digit stripped or suffixed tags (Fallback B).
- **Keys Stored in `global_mapping`**:
  - From `.gff3` files: `locus_tag` (or `ID`), `old_locus_tag`, and `gene`.
  - From `.xlsx` files: `locus tag`, `old locus tag`, `pan gene symbol`, `symbol`, `synonym`, and `queried_tag`.

---

## 5. Output and Artifact Generation

- **Dataset Modification**: The script inserts the new `gene` (symbol) column adjacent to the original locus tag column. `product`, `protein_id`, and `is_pseudogene` are appended to the dataset.
- **Summary Log**: A statistical summary (`add_gene_info_summary.xlsx`) is generated in the root directory, detailing the percentage of genes successfully annotated with products and protein IDs for every dataset.
- **FASTA Export**: If a gene is assigned a protein sequence but lacks a `WP_` RefSeq ID (mostly occurring with AureoWiki fallbacks), it is exported to `missing_protein_id.fasta`.
