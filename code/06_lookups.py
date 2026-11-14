"""Look up all instrument SNPs in the outcome GWAS and in the BMI / smoking-initiation GWAS.

  python3 -I code/06_lookups.py opengwas      OpenGWAS datasets, exact rsID matches (no proxies; cached per dataset)
  python3 -I code/06_lookups.py proxies       OpenGWAS proxy search (r2>=0.8) only for SNPs still missing
  python3 -I code/06_lookups.py repair        re-query SNPs still missing; rows the API returns under another rsID
                                              (dbSNP merges) are matched back to the queried rsID by GRCh37 chr:position
  python3 -I code/06_lookups.py local GCST..  stream a local GWAS Catalog harmonised file and keep the needed rsIDs
  python3 -I code/06_lookups.py expos         PSA_gw instruments in the exposure GWAS of replicated claims (MVMR; exact)

Extra rsIDs (e.g. PSA instruments) can be added through data/proc/extra_rsids.txt.
Outputs: work/lookups/<dataset>.tsv with columns rsid, ea, nea, eaf, beta, se, p, n, proxy, target_*.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
import gzip, time
import ogwas

OG_IDS = ["ieu-b-85", "ukb-b-19953", "ieu-b-4877", "ukb-b-15695", "ieu-b-4809", "finn-b-C3_PROSTATE_EXALLC"]
# ukb-b-19953 (UKB BMI, 9.8M variants) is used for the BMI correction instead of ieu-b-40 (HapMap3 only, 2.3M)
LDIR = os.path.join(WORK, "lookups")
os.makedirs(LDIR, exist_ok=True)
COLS = ["rsid", "chr", "position", "ea", "nea", "eaf", "beta", "se", "p", "n", "proxy", "target_snp", "proxy_snp", "target_a1", "target_a2", "proxy_a1", "proxy_a2"]


def needed_rsids():
    inst = pd.read_csv(os.path.join(PROC, "instruments.tsv"), sep="\t")
    rs = set(inst.rsid.astype(str))
    extra = os.path.join(PROC, "extra_rsids.txt")
    if os.path.exists(extra):
        rs |= {l.strip() for l in open(extra) if l.strip().startswith("rs")}
    return sorted(r for r in rs if r.startswith("rs"))


def load_cache(gid):
    p = os.path.join(LDIR, f"{gid}.tsv")
    if os.path.exists(p):
        return pd.read_csv(p, sep="\t", dtype={"rsid": str})
    return pd.DataFrame(columns=COLS + ["queried"])


def posmap():
    """GRCh37 chr:position of the exposure instruments (OpenGWAS coordinates)."""
    inst = pd.read_csv(os.path.join(PROC, "instruments.tsv"), sep="\t", dtype={"rsid": str, "chr": str})
    inst = inst.dropna(subset=["pos"]).drop_duplicates("rsid")
    return {r.rsid: (str(r.chr), int(r.pos)) for r in inst.itertuples()}


def opengwas(proxies=False, repair=False):
    rs = needed_rsids()
    need = set(rs)
    pm = posmap()
    for gid in OG_IDS:
        cache = load_cache(gid)
        if len(cache):
            cache = cache[cache.rsid.isin(need)]  # rows returned under another rsID are re-matched below
        done = set(cache.loc[cache.get("queried", pd.Series(dtype=str)).astype(str) == ("proxy" if proxies else "exact"), "rsid"]) if len(cache) else set()
        found = set(cache.loc[cache.beta.notna(), "rsid"]) if len(cache) else set()
        if repair:
            todo = [r for r in rs if r not in found and r in pm]
        elif proxies:
            todo = [r for r in rs if r not in found and r not in done]
        else:
            todo = [r for r in rs if r not in done]
        print(gid, "to query:", len(todo), "proxies" if proxies else ("repair" if repair else "exact"), flush=True)
        step = 64 if not proxies else 32
        new = []
        for i in range(0, len(todo), step):
            chunk = todo[i:i + step]
            cset = set(chunk)
            inv = {pm[r]: r for r in chunk if r in pm}
            data = {"id": [gid], "variant": chunk, "proxies": 1 if proxies else 0}
            if proxies:
                data.update({"r2": 0.8, "align_alleles": 1, "palindromes": 1, "maf_threshold": 0.3})
            res = ogwas.call("associations", data)
            got = set()
            for x in res:
                tgt = x.get("target_snp") or x.get("rsid")
                if tgt not in cset:
                    try:
                        tgt = inv.get((str(x.get("chr")), int(x.get("position"))))
                    except (TypeError, ValueError):
                        tgt = None
                    if tgt is None:
                        continue
                row = {c: x.get(c) for c in COLS}
                row["rsid"] = tgt
                row["returned_rsid"] = x.get("rsid")
                row["queried"] = "proxy" if proxies else "exact"
                new.append(row)
                got.add(tgt)
            for r in chunk:
                if r not in got:
                    new.append({"rsid": r, "queried": "proxy" if proxies else "exact"})
            if (i // step) % 20 == 0:
                print(f"  {gid} {i + len(chunk)}/{len(todo)} allowance={ogwas.LAST.get('remaining')}", flush=True)
                pd.concat([cache, pd.DataFrame(new)], ignore_index=True).to_csv(os.path.join(LDIR, f"{gid}.tsv"), sep="\t", index=False)
        cache = pd.concat([cache, pd.DataFrame(new)], ignore_index=True)
        # keep found rows preferentially
        cache["has"] = cache.beta.notna()
        cache = cache.sort_values(["rsid", "has"], ascending=[True, False]).drop_duplicates("rsid", keep="first").drop(columns="has")
        cache.to_csv(os.path.join(LDIR, f"{gid}.tsv"), sep="\t", index=False)
        print(gid, "found", int(cache.beta.notna().sum()), "of", len(rs), flush=True)


def local(gcst):
    path = os.path.join(RAW, "gwas", f"{gcst}.h.tsv.gz")
    rs = set(needed_rsids())
    out = []
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        irs = ix["rsid"]
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[irs] in rs:
                out.append(p)
    d = pd.DataFrame(out, columns=hdr)
    d = d.rename(columns={"effect_allele": "ea", "other_allele": "nea", "effect_allele_frequency": "eaf", "standard_error": "se",
                          "p_value": "p", "chromosome": "chr38", "base_pair_location": "pos38"})
    for c in ("beta", "se", "eaf", "p"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    if "beta" not in d or d.beta.isna().all():
        raise SystemExit("no beta column")
    # multi-allelic duplicates: keep the biallelic SNV record when possible
    d["snv"] = (d.ea.str.len() == 1) & (d.nea.str.len() == 1)
    d = d.sort_values(["rsid", "snv", "p"], ascending=[True, False, True]).drop_duplicates("rsid")
    d["proxy"] = False
    d[["rsid", "ea", "nea", "eaf", "beta", "se", "p", "proxy", "chr38", "pos38"]].to_csv(os.path.join(LDIR, f"{gcst}.tsv"), sep="\t", index=False)
    print(gcst, "found", len(d), "of", len(rs))


def expos():
    psa = pd.read_csv(os.path.join(PROC, "psa_instruments.tsv"), sep="\t", dtype={"rsid": str})
    rs = sorted(set(psa.loc[psa.key == "PSA_gw", "rsid"]))
    cm = pd.read_csv(os.path.join(PROC, "claims_map.tsv"), sep="\t")
    # MVMR is run for the claims that replicated in the primary outcome (results/tables/regrade_PSA_KLK.tsv)
    rg = pd.read_csv(os.path.join(TAB, "regrade_PSA_KLK.tsv"), sep="\t")
    keep = set(rg.loc[rg.replicated, "claim"])
    ids = sorted(set(cm.loc[(cm.inst_mode == "tophits") & cm.claim.isin(keep), "exp_id"]))
    edir = os.path.join(WORK, "lookups_exp")
    os.makedirs(edir, exist_ok=True)
    for gid in ids:
        p = os.path.join(edir, f"{gid}.tsv")
        if os.path.exists(p):
            continue
        res = ogwas.associations([gid], rs, proxies=0)
        d = pd.DataFrame([{c: x.get(c) for c in ("rsid", "ea", "nea", "eaf", "beta", "se", "p", "n")} for x in res])
        d.to_csv(p, sep="\t", index=False)
        print(gid, "found", len(d), "of", len(rs), "allowance", ogwas.LAST.get("remaining"), flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "opengwas":
        opengwas(False)
    elif mode == "proxies":
        opengwas(True)
    elif mode == "repair":
        opengwas(False, repair=True)
    elif mode == "local":
        local(sys.argv[2])
    elif mode == "expos":
        expos()
