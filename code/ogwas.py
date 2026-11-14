"""Minimal OpenGWAS API client. The token is read from .opengwas_token in the project root and is never printed or saved elsewhere."""
import os, json, time, urllib.request, urllib.parse
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
API = "https://api.opengwas.io/api"
_TOK = None
LAST = {}


def _token():
    global _TOK
    if _TOK is None:
        with open(os.path.join(_ROOT, ".opengwas_token")) as f:
            _TOK = f.read().strip()
    return _TOK


def call(endpoint, data=None, method=None, tries=5):
    url = f"{API}/{endpoint}"
    body = None
    if data is not None:
        body = urllib.parse.urlencode(data, doseq=True).encode()
    for i in range(tries):
        req = urllib.request.Request(url, data=body, method=method or ("POST" if body else "GET"),
                                     headers={"Authorization": "Bearer " + _token()})
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                LAST["remaining"] = r.headers.get("X-Allowance-Remaining")
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 429:  # allowance used up: wait for the reset announced by the server (at most 1 h), not a retry
                wait = min(3600, int(e.headers.get("Retry-After") or 600)) + 5
                print(f"OpenGWAS allowance exhausted; waiting {wait} s", flush=True)
                time.sleep(wait)
                return call(endpoint, data, method, tries)
            if e.code >= 500:
                time.sleep(10 * (i + 1)); continue
            raise RuntimeError(f"OpenGWAS {endpoint} HTTP {e.code}: {e.read()[:300]!r}")
        except Exception:
            if i == tries - 1: raise
            time.sleep(5 * (i + 1))
    raise RuntimeError(f"OpenGWAS {endpoint} failed after retries")


def gwasinfo(ids=None):
    if ids:
        return call("gwasinfo", {"id": list(ids)})
    return call("gwasinfo")


def tophits(gid, pval=5e-8, clump=1, r2=0.001, kb=10000):
    return call("tophits", {"id": gid, "pval": pval, "clump": clump, "r2": r2, "kb": kb, "preclumped": 0})


def associations(gids, rsids, proxies=0, r2=0.8, align_alleles=1, palindromes=1, maf_threshold=0.3):
    """The API allows N(id) x N(variant) <= 64 per request, so variants are chunked accordingly."""
    out = []
    gids, rsids = list(gids), list(rsids)
    step = max(1, 64 // len(gids))
    for i in range(0, len(rsids), step):
        d = {"id": gids, "variant": rsids[i:i + step], "proxies": proxies}
        if proxies:
            d.update({"r2": r2, "align_alleles": align_alleles, "palindromes": palindromes, "maf_threshold": maf_threshold})
        out += call("associations", d)
    return out
