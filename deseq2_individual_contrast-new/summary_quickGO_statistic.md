Ran locate QuickGO GAF file
Ran read top lines of GAF file
Ran python script to summarize GAF contents · Exited with code 1
Ran run robust summary of GAF · Exited with code 1
Ran check first 15 lines with utf-8-sig
Ran print columns of sample row
Ran analyze each column statistics
Ran check how GAF is used in ora_meta_analysis.R
Analyzed ora_meta_analysis.RL25-60
Ran check V2 and V3 samples

### Overview of the File

The file `QuickGO-annotations-1764241851574-20251127.gaf` is formatted in GAF 2.2 (Gene Association File format), exported from the EBI QuickGO  GOA (Gene Ontology Annotation) database on November 27, 2025.

#### Key Metadata & Dataset Summary
 Target Taxon Staphylococcus aureus (`taxon1280` and descendant strains such as N315, MW2, Mu50, and USA300).
 Total Records 804,442 annotation rows (and 11 comment lines starting with `!`).
 Unique GO Terms 2,344 distinct Gene Ontology IDs.
 Entities Annotated
   470,863 non-coding RNAs from [RNAcentral](httpsrnacentral.org) (58.5%)
   333,579 proteins from [UniProtKB](httpswww.uniprot.org) (41.5%)
 Ontology Aspect Distribution
   F (Molecular Function) 413,576 annotations (51.4%)
   P (Biological Process) 280,168 annotations (34.8%)
   C (Cellular Component) 110,698 annotations (13.8%)
 Evidence Quality
   IEA (Inferred from Electronic Annotation) 800,626 (99.5%)
   Curated  Phylogenetic  Experimental 3,816 annotations (including `IBA`, `ISS`, `IDA`, `IMP`, `EXP`, and `IGI`).

---

### Role of Each Column (GAF 2.2 Specification)

When imported into R using `read.delim(gaf_file, comment.char = !, header = FALSE)`, columns are indexed `V1` to `V17`

 Col (`R`)  Official GAF 2.2 Field  Required  Role & Description  Values in this Specific File 
 ---  ---  ---  ---  --- 
 `V1`  `DB`  Yes  Source Database Identifies the database repository managing the biological entity identifier.  `RNAcentral` (470,863) or `UniProtKB` (333,579). 
 `V2`  `DB_Object_ID`  Yes  Primary Accession  Identifier The unique accession code for the molecule in the source database.  UniProt accession (e.g., `A0A024BUM5`, `P0A0L1`) or RNAcentral ID (e.g., `URS000000007C_904785`). 
 `V3`  `DB_Object_Symbol`  Yes  Standard Gene  Molecule Symbol The commonly recognized symbol or label for the gene, protein, or RNA product.  Gene symbols (e.g., `mecA`, `pta`, `clfA`, `spa`, `agrA`) or full descriptive RNA labels (e.g., `Staphylococcus aureus ... tRNA-Lys`). 
 `V4`  `Qualifier`  Yes  Relation  Annotation Qualifier Formalized relation between the entity and the GO term (GAF 2.2 standard, with optional `NOT` prefix).  `enables` (MF terms 413,550), `involved_in` (BP terms 280,168), `located_in` (CC terms 102,705), `part_of`, `is_active_in`, `contributes_to`, `NOTenables`. 
 `V5`  `GO_ID`  Yes  Gene Ontology Identifier The standardized 7-digit accession ID for the attributed GO term.  e.g., `GO0008658`, `GO0005886`, `GO0071555`. 
 `V6`  `DBReference(s)`  Yes  Supporting Citation  Authority Publication or consortium reference supporting the annotation.  QuickGO reference IDs (`GO_REF0000120`, `GO_REF0000118`, `GO_REF0000002` [InterPro2GO], `GO_REF0000043` [KW2GO]) or PubMed IDs (`PMID...`). 
 `V7`  `Evidence_Code`  Yes  Evidence Code Standard 3-letter code indicating the experimental or computational origin of the assertion.  Predominantly `IEA` (automated electronic annotation 99.5%), followed by `IBA` (phylogenetic inference), `ISS` (sequence similarity), `IDA` (direct assay), and `EXP` (experimental). 
 `V8`  `With_From`  Conditional  Supporting Identifier Additional IDs referenced to justify electronic or similarity-based assertions.  InterPro domain signatures (e.g., `InterProIPR001460`), PANTHER families (`PANTHERPTN...`), or UniProt keywords (`UniProtKB-KW...`). 
 `V9`  `Aspect`  Yes  Ontology Namespace Single-letter code indicating which sub-ontology the term belongs to.  `F` = Molecular Functionbr`P` = Biological Processbr`C` = Cellular Component. 
 `V10`  `DB_Object_Name`  No  Full Gene  Protein Name Descriptive title of the molecule (e.g., Penicillin-binding protein 2a).  Empty (0%) in this QuickGO bulk download. 
 `V11`  `DB_Object_Synonym`  No  Synonyms  Aliases Pipe-separated alternative symbols or previous locus tags.  Empty (0%) in this QuickGO bulk download. 
 `V12`  `DB_Object_Type`  Yes  Entity Molecule Type Biological category of the annotated molecule.  `RNA` (470,863) or `protein` (333,579). 
 `V13`  `Taxon`  Yes  NCBI Taxonomy ID The organism NCBI Taxon ID (format `taxonID`).  `taxon1280` (S. aureus species 253,144) and specific substrains (e.g., `taxon93061` [N315], `taxon158878` [MW2], `taxon367830` [USA300]). 
 `V14`  `Date`  Yes  Date of Annotation Last updated  generated date in `YYYYMMDD` format.  Dates ranging from legacy entries up to recent curation releases (e.g., `20251107`). 
 `V15`  `Assigned_By`  Yes  Curating Authority Group or database pipeline asserting the annotation.  `RNAcentral`, `UniProt`, `GOC` (Gene Ontology Consortium), `TreeGrafter`, `InterPro`. 
 `V16`  `Annotation_Extension` No  Biological Context Extension Additional ontology clauses qualifying the annotation (substrates, locations, cell types).  Populated in 47 high-detail curated entries (e.g., `has_input(CHEBI4828)` or `occurs_in(...)`). 
 `V17`  `Gene_Product_Form_ID` No  Isoform  Modified Form ID Specific splice variant or PTM accession.  Empty (0%) — bacterial genomes do not utilize alternative splice variants. 

---

### Practical Considerations for Your Analysis (e.g., in [`ora_meta_analysis.R`](filecUsersAdminDocumentss_aueruss_auerus_biofilm_differential_projectdeseq2_individual_contrast-newora_meta_analysis.R))

In [`ora_meta_analysis.R#L30-L31`](filecUsersAdminDocumentss_aueruss_auerus_biofilm_differential_projectdeseq2_individual_contrast-newora_meta_analysis.R#L30-L31), you parse
```R
gaf - read.delim(gaf_file, comment.char = !, header = FALSE, stringsAsFactors = FALSE, quote = )
term2gene - unique(data.frame(term = gaf$V5, gene = gaf$V3, stringsAsFactors = FALSE))
```

Key takeaways to keep in mind
1. Filtering by Entity Type (`V12`) If your differential expression study focuses strictly on protein-coding genes, filter by `gaf$V12 == protein` (or `gaf$V1 == UniProtKB`) to prevent long RNAcentral labels (e.g., `Staphylococcus aureus... tRNA-Lys`) from diluting your gene mapping.
2. Filtering by Sub-ontology (`V9`) For pathwayprocess enrichment, you can split or filter by Biological Process (`gaf$V9 == P`), Molecular Function (`gaf$V9 == F`), or Cellular Component (`gaf$V9 == C`).
3. Symbol vs. Locus Tag Matching (`V2` vs `V3`) `gaf$V3` contains standard gene symbols (e.g., `mecA`), but uncharacterized proteins may default to UniProt accessions (e.g., `A0A059VJ41`). If your RNA-seq clusters use strain-specific locus tags (e.g., `SAUSA300_xxxx` or `SAVxxxx`), ensure your cluster mapping accounts for both gene symbols (`V3`) and UniProt accessions (`V2`).