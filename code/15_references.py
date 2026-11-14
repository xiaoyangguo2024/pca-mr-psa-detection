"""Reference list from PubMed records (NCBI E-utilities esummary), in AMA style (The Prostate) or Vancouver style.

Citation keys used in manuscript/manuscript.md as [@Key] or [@Key1; @Key2] map to PubMed IDs below. Records are cached
in lit/pubmed_esummary.json so the manuscript can be rebuilt offline. Only PMIDs are sent to NCBI.
Output: lit/pubmed_esummary.json
Run: python3 -I code/15_references.py
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
import re, time, urllib.request

KEYS = {
    "Cui2024": "38489391", "Bell2015": "25821151", "Banez2007": "18029831", "Discacciati2012": "22228452",
    "Dluzniewski2015": "25352524", "Kachuri2023": "37264206", "Shi2023": "37865044", "Hoffmann2025": "39930085",
    "Penney2011": "21421545", "Skrivankova2021": "34698778", "Hemani2018": "29846171", "Wang2023": "37945903",
    "Schumacher2018": "29892016", "Kurki2023": "36653562", "LloydJones2018": "29429966", "Aschard2015": "25640676",
    "Gilbody2025": "39812504", "WangP2025": "39891386", "Liu2019": "30643251", "Yavorska2017": "28398548",
    "Burgess2019": "32760811", "Bowden2015": "26050253", "Bowden2016": "27061298", "Hartwig2017": "29040600",
    "Islami2014": "25242554", "Watts2022": "35579976", "Burgess2016": "27625185", "Kazmi2020": "31802111",
}
PAGES = {"29846171": "e34408"}  # eLife article number (PubMed esummary gives no pages for this record)
CACHE = os.path.join(ROOT, "lit", "pubmed_esummary.json")


def fetch(pmids):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(pmids))
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)["result"]


def load():
    if os.path.exists(CACHE):
        c = load_json(CACHE)
        if all(p in c for p in KEYS.values()):
            return c
    ids, c = sorted(set(KEYS.values())), {}
    for i in range(0, len(ids), 4):  # small batches: long chunked responses are cut by the proxy
        for t in range(5):
            try:
                res = fetch(ids[i:i + 4])
                break
            except Exception as e:
                print("retry", t + 1, type(e).__name__)
                time.sleep(2)
        c.update({u: res[u] for u in res["uids"]})
        time.sleep(0.4)
    save_json(c, CACHE)
    return c


def vancouver(x):
    au = [a["name"] for a in x["authors"] if a.get("authtype") == "Author"]
    names = ", ".join(au[:6]) + (", et al." if len(au) > 6 else ".")
    title = x["title"].rstrip(".")
    year = x["pubdate"][:4]
    vol, iss, pg = x.get("volume", ""), x.get("issue", ""), x.get("pages", "")
    pg = pg or PAGES.get(x["uid"], "")
    if not pg:  # online-only articles: article number from the electronic location (e.g. "pii: e34408. doi: ...")
        m = re.search(r"\b(e\d+)\b", x.get("elocationid", ""))
        pg = m.group(1) if m else ""
    loc = f"{year};{vol}" + (f"({iss})" if iss else "") + (f":{pg}" if pg else "")
    doi = next((a["value"] for a in x.get("articleids", []) if a["idtype"] == "doi"), "")
    return f"{names} {title}. {x['source']}. {loc}." + (f" doi:{doi}" if doi else "")


def full_pages(pg):
    """Expand abbreviated page ranges (1749-57 -> 1749-1757)."""
    m = re.fullmatch(r"(\d+)-(\d+)", pg)
    if not m or len(m.group(2)) >= len(m.group(1)):
        return pg
    a, b = m.groups()
    return f"{a}-{a[:len(a) - len(b)]}{b}"


def ama(x):
    """AMA style as in The Prostate's sample: 'Author AB, Author CD. Title. J Abbrev Year;Vol:pages.'
    All authors when six or fewer, otherwise the first three followed by et al."""
    au = [a["name"] for a in x["authors"] if a.get("authtype") == "Author"]
    names = ", ".join(au) if len(au) <= 6 else ", ".join(au[:3]) + ", et al"
    title = x["title"].rstrip(".")
    src = {"Elife": "eLife"}.get(x["source"], x["source"])
    pg = full_pages(x.get("pages", "") or PAGES.get(x["uid"], ""))
    return f"{names}. {title}. {src} {x['pubdate'][:4]};{x.get('volume', '')}" + (f":{pg}" if pg else "") + "."


if __name__ == "__main__":
    c = load()
    for k, p in KEYS.items():
        print(k, "|", ama(c[p]))
