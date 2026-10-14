import requests
from bs4 import BeautifulSoup
import re

def search_gene_api(gene_symbol):
    """
    Uses the MediaWiki API to search for a gene and return the top matching page titles.
    """
    api_url = "https://aureowiki.med.uni-greifswald.de/api.php"
    
    params = {
        "action": "query",
        "list": "search",
        "srsearch": gene_symbol,
        "utf8": "1",
        "format": "json",
        "srlimit": 20  # Fetch the top 20 results
    }
    
    response = requests.get(api_url, params=params, headers={"User-Agent": "AureoWikiQueryBot/1.0"})
    data = response.json()
    
    search_results = data.get("query", {}).get("search", [])
    
    if not search_results:
        print(f"[-] No results found for '{gene_symbol}'.")
        return None
        
    print(f"[+] Found {len(search_results)} top results for '{gene_symbol}':")
    
    # Store just the exact page titles
    page_titles = []
    for idx, result in enumerate(search_results):
        title = result["title"]
        page_titles.append(title)
        print(f"    {idx + 1}. {title}")
        
    # Return the first matching title (or you can modify this to filter by a specific strain)
    return page_titles[0]


def fetch_rendered_data(page_title):
    """
    Fetches the rendered HTML for the specific locus tag/page title.
    """
    api_url = "https://aureowiki.med.uni-greifswald.de/api.php"
    params = {
        "action": "parse",
        "page": page_title,
        "prop": "text",
        "format": "json"
    }
    
    response = requests.get(api_url, params=params, headers={"User-Agent": "AureoWikiQueryBot/1.0"})
    data = response.json()
    
    if "parse" in data and "text" in data["parse"]:
        html_content = data["parse"]["text"]["*"]
        soup = BeautifulSoup(html_content, "html.parser")
        
        print(f"\n--- Data for {page_title} ---")
        summary_header = soup.find(lambda tag: tag.name in ["h1", "h2"] and "Summary" in tag.text)
        if summary_header:
            ul = summary_header.find_next("ul")
            if ul:
                for li in ul.find_all("li", recursive=False):
                    print(li.get_text(strip=True))
    else:
        print("Failed to fetch rendered data.")

def fetch_protein_accession(page_title):
    """
    Fetches the rendered HTML and extracts the RefSeq protein accession number (e.g., YP_, WP_).
    """
    api_url = "https://aureowiki.med.uni-greifswald.de/api.php"
    params = {
        "action": "parse",
        "page": page_title,
        "prop": "text",
        "format": "json"
    }
    
    response = requests.get(api_url, params=params, headers={"User-Agent": "AureoWikiQueryBot/1.0"})
    data = response.json()
    
    if "parse" in data and "text" in data["parse"]:
        html_content = data["parse"]["text"]["*"]
        soup = BeautifulSoup(html_content, "html.parser")
        
        print(f"\n--- Extracting Protein Accession for {page_title} ---")
        
        # Use a separator to ensure spaces between nested tags (like <a> links) are preserved
        for li in soup.find_all("li"):
            text = li.get_text(separator=" ", strip=True)
            
            # The Accession numbers list item typically looks like "RefSeq: YP_494378 NCBI"
            if text.startswith("RefSeq:"):
                # Use regex to find typical NCBI protein prefixes (YP_, WP_, NP_, MGT, etc.)
                match = re.search(r'(YP_\d+|WP_\d+|NP_\d+|MGT\d+)', text)
                
                if match:
                    protein_accession = match.group(1)
                    print(f"Raw list item: {text}")
                    print(f"[+] Extracted Protein Accession: {protein_accession}")
                    return protein_accession
                    
        print("[-] No RefSeq protein accession found on this page.")
        
    else:
        print("Failed to fetch rendered data.")

# --- Run the test ---
gene_to_search = "purH"
top_result_title = search_gene_api(gene_to_search)

if top_result_title:
    print(f"\n[*] Automatically fetching data for the top result: {top_result_title}")
    fetch_rendered_data(top_result_title)
    fetch_protein_accession(top_result_title)