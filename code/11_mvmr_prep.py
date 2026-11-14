"""Inputs for multivariable MR: exposure X and ln PSA (BMI/smoking adjustment undone) on prostate cancer (Wang 2023).

For each genome-wide claim: SNP set = X instruments + genome-wide PSA instruments, clumped together with plink on the
1000 Genomes EUR reference (OpenGWAS ld/clump; r2 < 0.001, 10 Mb; index p = the SNP's p-value for the exposure it was
selected for).
Effects are aligned to the allele of the list the SNP came from:
  X SNPs:   bX from the X GWAS (instrument row), bP and bY from data/proc/harmonised.tsv (same key)
  PSA SNPs: bP (corrected, outcome PSA_total) and bY from harmonised.tsv (key PSA_gw), bX from work/lookups_exp/<exp_id>.tsv (exact rsID match)
Run for the genome-wide claims that replicated in the primary outcome; drug-target (region) claims are not analysed.
Output: data/proc/mvmr_input.tsv
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
import ogwas

h = pd.read_csv(os.path.join(PROC, "harmonised.tsv"), sep="\t", dtype={"rsid": str})
inst = pd.read_csv(os.path.join(PROC, "instruments.tsv"), sep="\t", dtype={"rsid": str})
psa = pd.read_csv(os.path.join(PROC, "psa_instruments.tsv"), sep="\t", dtype={"rsid": str})
psa = psa[psa.key == "PSA_gw"]
cm = pd.read_csv(os.path.join(PROC, "claims_map.tsv"), sep="\t")


def wide(key, rs):
    s = h[(h.key == key) & h.rsid.isin(rs)]
    p = s[s.outcome == "PSA_total"].set_index("rsid")[["ea", "nea", "eaf", "bx", "sx", "by", "sy"]].rename(columns={"by": "bP", "sy": "sP"})
    y = s[s.outcome == "PCa_Wang2023"].set_index("rsid")[["by", "sy"]].rename(columns={"by": "bY", "sy": "sY"})
    return p.join(y, how="inner")


def align(ref, look):
    """ref: rows with ea, nea, eaf; look: GWAS rows (rsid, ea, nea, eaf, beta, se). Returns beta, se aligned to ref.ea."""
    m = ref.join(look.set_index("rsid"), rsuffix="_l", how="left")
    out_b, out_s = [], []
    for r in m.itertuples():
        ea, nea, la, lna = str(r.ea).upper(), str(r.nea).upper(), str(r.ea_l).upper(), str(r.nea_l).upper()
        if pd.isna(r.beta):
            out_b.append(np.nan); out_s.append(np.nan); continue
        if is_palindromic(ea, nea):
            if pd.isna(r.eaf) or pd.isna(r.eaf_l) or min(r.eaf, 1 - r.eaf) > 0.42:
                out_b.append(np.nan); out_s.append(np.nan); continue
            sg = 1.0 if (r.eaf - 0.5) * (r.eaf_l - 0.5) > 0 else -1.0
        elif (la, lna) in ((ea, nea), (COMP.get(ea), COMP.get(nea))):
            sg = 1.0
        elif (la, lna) in ((nea, ea), (COMP.get(nea), COMP.get(ea))):
            sg = -1.0
        else:
            out_b.append(np.nan); out_s.append(np.nan); continue
        out_b.append(sg * r.beta); out_s.append(r.se)
    return np.array(out_b, float), np.array(out_s, float)


rows = []
rg = pd.read_csv(os.path.join(TAB, "regrade_PSA_KLK.tsv"), sep="\t")
for c in cm[(cm.inst_mode == "tophits") & cm.claim.isin(set(rg.loc[rg.replicated, "claim"]))].itertuples():
    key = f"claim{c.claim:02d}"
    xi = inst[inst.key == key]
    pool = pd.concat([xi[["rsid", "p"]].assign(src="X"), psa[["rsid", "p"]].assign(src="PSA")]).sort_values("p")
    pool = pool.drop_duplicates("rsid")
    keep = set(ogwas.call("ld/clump", {"rsid": list(pool.rsid), "pval": [float(x) for x in pool.p], "pthresh": 1.0,
                                       "r2": 0.001, "kb": 10000, "pop": "EUR"}))
    pool = pool[pool.rsid.isin(keep)]
    xs = wide(key, set(pool.loc[pool.src == "X", "rsid"]))
    xs["bX"], xs["sX"] = xs.bx, xs.sx
    ps = wide("PSA_gw", set(pool.loc[pool.src == "PSA", "rsid"]))
    lp = os.path.join(WORK, "lookups_exp", f"{c.exp_id}.tsv")
    if not os.path.exists(lp) or ps.empty:
        continue
    look = pd.read_csv(lp, sep="\t", dtype={"rsid": str})
    for col in ("beta", "se", "eaf"):
        look[col] = pd.to_numeric(look[col], errors="coerce")
    ps["bX"], ps["sX"] = align(ps, look)
    both = pd.concat([xs.assign(src="X"), ps.assign(src="PSA")])
    both = both[np.isfinite(both[["bX", "sX", "bP", "sP", "bY", "sY"]]).all(axis=1)]
    both["key"] = key
    rows.append(both.reset_index()[["key", "rsid", "src", "bX", "sX", "bP", "sP", "bY", "sY"]])
    print(key, "X SNPs", int((both.src == "X").sum()), "PSA SNPs", int((both.src == "PSA").sum()), flush=True)
pd.concat(rows).to_csv(os.path.join(PROC, "mvmr_input.tsv"), sep="\t", index=False)
