"""Small GWAS Catalog REST v2 helpers (studies by PubMed ID, curated associations of a study)."""
import json, time, urllib.request, urllib.parse
API = "https://www.ebi.ac.uk/gwas/rest/api/v2"

def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            if i == tries - 1: raise
            time.sleep(2 + 3 * i)

def studies_by_pmid(pmid):
    out, page = [], 0
    while True:
        d = get(f"{API}/studies?pubmed_id={pmid}&size=200&page={page}")
        out += d.get("_embedded", {}).get("studies", [])
        if page + 1 >= d.get("page", {}).get("totalPages", 1): break
        page += 1
    return out

def associations(acc):
    out, page = [], 0
    while True:
        d = get(f"{API}/associations?accession_id={acc}&size=500&page={page}")
        out += d.get("_embedded", {}).get("associations", [])
        if page + 1 >= d.get("page", {}).get("totalPages", 1): break
        page += 1
    return out
