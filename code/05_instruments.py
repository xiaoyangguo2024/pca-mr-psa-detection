"""Instrument extraction from OpenGWAS for every re-analysable claim, plus BMI and smoking initiation (used to undo the
covariate adjustment of the PSA GWAS).

tophits mode: server-side clumping (EUR 1000G, r2<0.001, 10 Mb) at the claim's threshold; when the threshold is 5e-8 and
fewer than 3 SNPs remain, it is relaxed to 5e-6 and then 1e-5 (the threshold used is recorded).
region mode (drug targets): all variants within +-100 kb of the target genes (GRCh37 Ensembl coordinates) with p<5e-8,
clumped at r2<0.1 within 250 kb (server-side), plus the EUR LD matrix for a correlation-aware IVW.
Outputs (cached, re-run safe): work/instruments/<key>.json, data/proc/instruments.tsv, data/proc/instrument_summary.tsv
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
import ogwas

GENES37 = {  # Ensembl GRCh37 lookup, 2026-11-13
    "PCSK9": ("1", 55505221, 55530525), "PPARG": ("3", 12328867, 12475855),
    "CACNA1C": ("12", 2079952, 2802108), "CACNA1D": ("3", 53528683, 53847760), "CACNA1F": ("X", 49061523, 49089833),
    "CACNA1S": ("1", 201008642, 201081694), "CACNA2D1": ("7", 81575760, 82073114), "CACNA2D2": ("3", 50400233, 50541675),
    "CACNB1": ("17", 37329709, 37353956), "CACNB2": ("10", 18429606, 18830798), "CACNB3": ("12", 49207577, 49222724),
    "CACNB4": ("2", 152689290, 152955593), "CACNG1": ("17", 65040706, 65052909),
}
IDIR = os.path.join(WORK, "instruments")
os.makedirs(IDIR, exist_ok=True)


def cached(key, fn):
    p = os.path.join(IDIR, key + ".json")
    if os.path.exists(p):
        return load_json(p)
    obj = fn()
    save_json(obj, p)
    return obj


def tophits_relaxed(gid, p0):
    used, hits = p0, []
    for p in ([p0] if p0 != 5e-8 else [5e-8, 5e-6, 1e-5]):
        hits = ogwas.tophits(gid, pval=p)
        used = p
        if len(hits) >= 3:
            break
    return dict(p_used=used, hits=hits)


def region_hits(gid, genes):
    rows = []
    for g in genes:
        ch, s, e = GENES37[g]
        rng = f"{ch}:{max(1, s - 100000)}-{e + 100000}"
        a = ogwas.call("associations", {"id": [gid], "variant": [rng], "proxies": 0})
        for x in a:
            x["gene"] = g
        rows += a
    rows = [x for x in rows if x.get("p") is not None and float(x["p"]) < 5e-8 and str(x.get("rsid", "")).startswith("rs")]
    if not rows:
        return dict(p_used=5e-8, hits=[], ld=None)
    # keep one record per rsid
    seen = {}
    for x in rows:
        if x["rsid"] not in seen or float(x["p"]) < float(seen[x["rsid"]]["p"]):
            seen[x["rsid"]] = x
    rows = list(seen.values())
    cl = ogwas.call("ld/clump", {"rsid": [x["rsid"] for x in rows], "pval": [float(x["p"]) for x in rows],
                                 "pthresh": 5e-8, "r2": 0.1, "kb": 250, "pop": "EUR"})
    keep = set(cl)
    hits = [x for x in rows if x["rsid"] in keep]
    ld = ogwas.call("ld/matrix", {"rsid": [x["rsid"] for x in hits], "pop": "EUR"}) if len(hits) > 1 else None
    return dict(p_used=5e-8, hits=hits, ld=ld)


cm = pd.read_csv(os.path.join(PROC, "claims_map.tsv"), sep="\t")
jobs = [("BMI", "ieu-b-40", "tophits", 5e-8, ""), ("SMK", "ieu-b-4877", "tophits", 5e-8, "")]
for r in cm.itertuples():
    if r.inst_mode == "exclude":
        continue
    jobs.append((f"claim{r.claim:02d}", r.exp_id, r.inst_mode, float(r.inst_p), r.region_genes if isinstance(r.region_genes, str) else ""))

rows, summ = [], []
for key, gid, mode, p0, genes in jobs:
    if mode == "tophits":
        obj = cached(key, lambda: tophits_relaxed(gid, p0))
    else:
        obj = cached(key, lambda: region_hits(gid, genes.split(",")))
    hits = obj["hits"]
    for h in hits:
        rows.append(dict(key=key, exp_id=gid, rsid=h["rsid"], chr=str(h.get("chr")), pos=h.get("position"), ea=h.get("ea"),
                         nea=h.get("nea"), eaf=h.get("eaf"), beta=h.get("beta"), se=h.get("se"), p=h.get("p"), n=h.get("n"),
                         gene=h.get("gene", "")))
    summ.append(dict(key=key, exp_id=gid, mode=mode, p_used=obj["p_used"], n_snp=len(hits)))
    print(key, gid, mode, "p", obj["p_used"], "n", len(hits), "| allowance left", ogwas.LAST.get("remaining"), flush=True)

inst = pd.DataFrame(rows)
for c in ("beta", "se", "p", "eaf"):
    inst[c] = pd.to_numeric(inst[c], errors="coerce")
inst["F"] = (inst.beta / inst.se) ** 2
inst.to_csv(os.path.join(PROC, "instruments.tsv"), sep="\t", index=False)
s = pd.DataFrame(summ)
fs = inst.groupby("key").F.agg(["mean", "min"]).rename(columns={"mean": "F_mean", "min": "F_min"})
s = s.merge(fs, left_on="key", right_index=True, how="left")
s.to_csv(os.path.join(PROC, "instrument_summary.tsv"), sep="\t", index=False)
print(s.to_string(index=False))
