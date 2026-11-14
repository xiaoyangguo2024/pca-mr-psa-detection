"""PSA-level instruments (for theta, the effect of genetically higher PSA on prostate cancer detection).

Source: Hoffmann 2025 PSA GWAS (GCST90461907, harmonised GRCh38). Rows with p < 5e-8, biallelic SNVs with an rsID.
  PSA_gw   genome-wide: clumping on the 1000 Genomes EUR reference (OpenGWAS ld/clump, plink-based; r2 < 0.001, 10 Mb),
           run per chromosome
  PSA_KLK  cis-KLK2/KLK3 locus (GRCh38 chr19:50.70-51.00 Mb): clumping r2 < 0.1 within 250 kb, plus the EUR LD matrix
           (OpenGWAS ld/matrix, EUR, signed for the first listed allele) stored in work/instruments/PSA_KLK.json for
           correlation-aware IVW.
           KLK3 encodes PSA itself, so these variants shift measured PSA with the least scope for other biology.
Outputs: data/proc/psa_instruments.tsv, data/proc/extra_rsids.txt (appended), work/instruments/PSA_KLK.json
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
import gzip
import ogwas

SRC = os.path.join(RAW, "gwas", "GCST90461907.h.tsv.gz")
KLK = ("19", 50_700_000, 51_000_000)  # GRCh38

rows = []
with gzip.open(SRC, "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {h: i for i, h in enumerate(hdr)}
    ip = ix["p_value"]
    for line in f:
        p = line.rstrip("\n").split("\t")
        try:
            pv = float(p[ip])
        except ValueError:
            continue
        if pv < 5e-8:
            rows.append(p)
d = pd.DataFrame(rows, columns=hdr)
print("columns:", hdr)
print("p<5e-8 rows:", len(d))
d = d.rename(columns={"effect_allele": "ea", "other_allele": "nea", "effect_allele_frequency": "eaf", "standard_error": "se",
                      "p_value": "p", "chromosome": "chr", "base_pair_location": "pos"})
for c in ("beta", "se", "eaf", "p"):
    d[c] = pd.to_numeric(d[c], errors="coerce")
d["pos"] = pd.to_numeric(d["pos"], errors="coerce")
d = d[d.rsid.astype(str).str.startswith("rs") & (d.ea.str.len() == 1) & (d.nea.str.len() == 1) & d.beta.notna()]
d = d.sort_values("p").drop_duplicates("rsid")
ncol = "n" if "n" in d.columns else None
print("biallelic SNVs with rsID:", len(d))


def clump(sub, tag, r2, kb):
    keep = []
    for ch, g in sub.groupby("chr"):
        res = ogwas.call("ld/clump", {"rsid": list(g.rsid), "pval": [float(x) for x in g.p], "pthresh": 5e-8, "r2": r2,
                                      "kb": kb, "pop": "EUR"})
        keep += list(res)
        print(f"  {tag} chr{ch}: {len(g)} -> {len(res)}", flush=True)
    return keep


gw = clump(d, "psa_gw", 0.001, 10000)
klk = d[(d.chr == KLK[0]) & d.pos.between(KLK[1], KLK[2])]
kl = clump(klk, "psa_klk", 0.1, 250) if len(klk) else []
print("PSA_gw instruments:", len(gw), "| PSA_KLK instruments:", len(kl), "of", len(klk), "KLK rows")

# LD matrix for the KLK set
ld = ogwas.call("ld/matrix", {"rsid": kl, "pop": "EUR"}) if len(kl) > 1 else None
save_json(dict(p_used=5e-8, hits=[], ld=ld, note="cis-KLK PSA instruments; LD from OpenGWAS ld/matrix (1000G EUR)"),
          os.path.join(WORK, "instruments", "PSA_KLK.json"))

out = []
for key, ids in (("PSA_gw", gw), ("PSA_KLK", kl)):
    s = d[d.rsid.isin(ids)].copy()
    s["key"] = key
    s["exp_id"] = "GCST90461907"
    s["n"] = pd.to_numeric(s[ncol], errors="coerce") if ncol else np.nan
    s["gene"] = "KLK" if key == "PSA_KLK" else ""
    out.append(s[["key", "exp_id", "rsid", "chr", "pos", "ea", "nea", "eaf", "beta", "se", "p", "n", "gene"]])
psa = pd.concat(out, ignore_index=True)
psa["F"] = (psa.beta / psa.se) ** 2
psa.to_csv(os.path.join(PROC, "psa_instruments.tsv"), sep="\t", index=False)
extra = os.path.join(PROC, "extra_rsids.txt")
old = set(open(extra).read().split()) if os.path.exists(extra) else set()
open(extra, "w").write("\n".join(sorted(old | set(psa.rsid))) + "\n")
print(psa.groupby("key").agg(n=("rsid", "size"), F_mean=("F", "mean"), F_min=("F", "min")))
