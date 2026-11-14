"""Harmonise exposure instruments with every outcome and build the BMI/smoking-corrected PSA outcome.

Exposure sets: data/proc/instruments.tsv (claims, BMI, SMK) and data/proc/psa_instruments.tsv (PSA level, for theta).
Outcome lookups: work/lookups/<dataset>.tsv (06_lookups.py).

Allele rules: match or flip on the forward strand, then on the complementary strand; palindromic SNPs are kept only when
both allele frequencies are available and the exposure MAF <= 0.42, and are oriented by frequency; others dropped.
OpenGWAS proxy rows are translated to the target SNP's alleles with target_a1/a2 <-> proxy_a1/a2.

PSA_total (the PSA GWAS adjusted for BMI and smoking undone, Aschard 2015 AJHG):
    b_total = b_adj + c_BMI * SD_BMI * b_G,BMI(SD) + c_SMK * b_G,SMK(log-odds) * p(1-p)
c_BMI from NHANES (work/nhanes_c.json, per kg/m2); SD_BMI = 4.8 kg/m2 (UK Biobank, the BMI lookup GWAS ukb-b-19953);
smoking: ever-smoking coefficient applied to the change in prevalence implied by the smoking-initiation log-odds (p=0.5).
Output: data/proc/harmonised.tsv (one row per exposure set x outcome x SNP) and data/proc/harmonise_log.tsv
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *

LDIR = os.path.join(WORK, "lookups")
OUTCOMES = {  # outcome label -> lookup file
    "PCa_Wang2023": "GCST90274714", "PCa_PRACTICAL": "ieu-b-85", "PCa_UKB": "ieu-b-4809",
    "PCa_FinnGen": "finn-b-C3_PROSTATE_EXALLC", "PSA_adj": "GCST90461907", "PSAtest_UKB": "ukb-b-15695",
    "BMI_UKB": "ukb-b-19953", "SMK_init": "ieu-b-4877",
}
# GWAS run as linear models on 0/1 outcomes (BOLT-LMM): beta and se are on the risk-difference scale and are converted
# to the log-odds scale with beta / (mu * (1 - mu)), mu = case fraction (checked: the slope of ieu-b-4809 SNP effects on
# Wang 2023 log ORs for 57 strong SNPs was 0.0486 vs mu(1-mu) = 0.0475).
LINEAR_OUT = {"ieu-b-4809": 9132 / 182625, "ukb-b-15695": 62489 / 200410}
LINEAR_EXP = {"ukb-b-20289": 3545 / 462933}
SD_BMI = 4.8
P_SMK = 0.5
NH = load_json(os.path.join(WORK, "nhanes_c.json"))
NHT = pd.read_csv(os.path.join(TAB, "nhanes_psa_covariates.tsv"), sep="\t")
C_EVER = float(NHT[(NHT.model == "BMI + ever smoking") & (NHT.term == "ever_smk")].beta.iloc[0])


def comp(s):
    return s.map(lambda a: "".join(COMP.get(c, "?") for c in a) if isinstance(a, str) else a)


def load_outcome(fid):
    d = pd.read_csv(os.path.join(LDIR, f"{fid}.tsv"), sep="\t", dtype={"rsid": str, "ea": str, "nea": str})
    d = d[d.beta.notna()].copy()
    for c in ("ea", "nea"):
        d[c] = d[c].str.upper()
    if "target_a1" in d.columns:
        px = d.proxy.astype(str).str.lower().isin(["true", "1"])
        for c in ("target_a1", "target_a2", "proxy_a1", "proxy_a2"):
            d[c] = d[c].astype(str).str.upper()
        # translate proxy alleles to target alleles
        ea_t = np.where(d.ea == d.proxy_a1, d.target_a1, np.where(d.ea == d.proxy_a2, d.target_a2, None))
        nea_t = np.where(d.nea == d.proxy_a1, d.target_a1, np.where(d.nea == d.proxy_a2, d.target_a2, None))
        d.loc[px, "ea"] = ea_t[px]
        d.loc[px, "nea"] = nea_t[px]
        d = d[~(px & (pd.isna(d.ea) | pd.isna(d.nea)))]
        d["proxy"] = px
    else:
        d["proxy"] = False
    if fid in LINEAR_OUT:
        k = LINEAR_OUT[fid] * (1 - LINEAR_OUT[fid])
        d["beta"], d["se"] = d.beta / k, d.se / k
    keep = ["rsid", "ea", "nea", "eaf", "beta", "se", "p", "proxy"]
    return d[keep].drop_duplicates("rsid")


def harmonise(x, y):
    """x: exposure rows (rsid, ea, nea, eaf, beta, se); y: outcome rows. Returns merged rows with by aligned to x.ea."""
    m = x.merge(y, on="rsid", suffixes=("_x", "_y"))
    if m.empty:
        return m.assign(action=[], by=[], palin=[])
    xa, xb = m.ea_x.str.upper(), m.nea_x.str.upper()
    ya, yb = m.ea_y, m.nea_y
    same = (ya == xa) & (yb == xb)
    swap = (ya == xb) & (yb == xa)
    csame = (comp(ya) == xa) & (comp(yb) == xb)
    cswap = (comp(ya) == xb) & (comp(yb) == xa)
    sign = np.where(same | csame, 1.0, np.where(swap | cswap, -1.0, np.nan))
    palin = pd.Series([is_palindromic(a, b) for a, b in zip(xa, xb)], index=m.index)
    exf = pd.to_numeric(m.eaf_x, errors="coerce")
    eyf = pd.to_numeric(m.eaf_y, errors="coerce")
    maf = np.minimum(exf, 1 - exf)
    ok_pal = palin & exf.notna() & eyf.notna() & (maf <= 0.42)
    # palindromic: orient by frequency (y eaf refers to y.ea; same side of 0.5 -> same allele)
    pal_sign = np.where((exf - 0.5) * (eyf - 0.5) > 0, 1.0, -1.0)
    # if y alleles were literally swapped relative to x, y.eaf refers to x.nea -> handled by pal rule directly
    sign = np.where(palin, np.where(ok_pal, pal_sign, np.nan), sign)
    m["action"] = np.where(np.isnan(sign), "drop_allele", np.where(palin, "palindrome_kept", "ok"))
    m.loc[palin & ~ok_pal, "action"] = "drop_palindrome"
    m["by"] = m.beta_y * sign
    m["palin"] = palin
    return m


def build(exp, label_col):
    rows, log = [], []
    outs = {k: load_outcome(v) for k, v in OUTCOMES.items() if os.path.exists(os.path.join(LDIR, f"{v}.tsv"))}
    for key, x in exp.groupby(label_col):
        x = x.drop_duplicates("rsid")
        per = {}
        for lab, y in outs.items():
            h = harmonise(x[["rsid", "ea", "nea", "eaf", "beta", "se", "p"]], y)
            log.append(dict(key=key, outcome=lab, n_inst=len(x), n_found=len(h), n_proxy=int(h.proxy.sum()) if len(h) else 0,
                            n_drop_allele=int((h.action == "drop_allele").sum()) if len(h) else 0,
                            n_drop_pal=int((h.action == "drop_palindrome").sum()) if len(h) else 0))
            h = h[h.action.isin(["ok", "palindrome_kept"])]
            per[lab] = h.set_index("rsid")
            for r in h.itertuples():
                rows.append(dict(key=key, outcome=lab, rsid=r.rsid, ea=r.ea_x, nea=r.nea_x, eaf=r.eaf_x, bx=r.beta_x, sx=r.se_x,
                                 px=r.p_x, by=r.by, sy=r.se_y, py=r.p_y, proxy=bool(r.proxy), palin=bool(r.palin)))
        # corrected PSA (total effect)
        if "PSA_adj" in per:
            pa = per["PSA_adj"]
            bmi = per.get("BMI_UKB", pd.DataFrame(columns=["by", "se_y"]))
            smk = per.get("SMK_init", pd.DataFrame(columns=["by", "se_y"]))
            for rs, r in pa.iterrows():
                bb = bmi.by.get(rs, np.nan) if len(bmi) else np.nan
                sb = bmi.se_y.get(rs, np.nan) if len(bmi) else np.nan
                bs = smk.by.get(rs, np.nan) if len(smk) else np.nan
                ss = smk.se_y.get(rs, np.nan) if len(smk) else np.nan
                kB = NH["c_bmi_per_kgm2"] * SD_BMI
                kS = C_EVER * P_SMK * (1 - P_SMK)
                tb = kB * (bb if pd.notna(bb) else 0.0) + kS * (bs if pd.notna(bs) else 0.0)
                tv = (kB * sb) ** 2 if pd.notna(sb) else 0.0
                tv += (kS * ss) ** 2 if pd.notna(ss) else 0.0
                rows.append(dict(key=key, outcome="PSA_total", rsid=rs, ea=r.ea_x, nea=r.nea_x, eaf=r.eaf_x, bx=r.beta_x,
                                 sx=r.se_x, px=r.p_x, by=r.by + tb, sy=math.sqrt(r.se_y ** 2 + tv), py=np.nan,
                                 proxy=bool(r.proxy), palin=bool(r.palin), bmi_term=kB * bb if pd.notna(bb) else np.nan,
                                 smk_term=kS * bs if pd.notna(bs) else np.nan))
    return pd.DataFrame(rows), pd.DataFrame(log)


if __name__ == "__main__":
    inst = pd.read_csv(os.path.join(PROC, "instruments.tsv"), sep="\t", dtype={"rsid": str})
    sets = [inst]
    pp = os.path.join(PROC, "psa_instruments.tsv")
    if os.path.exists(pp):
        sets.append(pd.read_csv(pp, sep="\t", dtype={"rsid": str}))
    allx = pd.concat(sets, ignore_index=True)
    for gid, mu in LINEAR_EXP.items():
        m = allx.exp_id == gid
        allx.loc[m, "beta"] = allx.loc[m, "beta"] / (mu * (1 - mu))
        allx.loc[m, "se"] = allx.loc[m, "se"] / (mu * (1 - mu))
    h, log = build(allx, "key")
    h.to_csv(os.path.join(PROC, "harmonised.tsv"), sep="\t", index=False)
    log.to_csv(os.path.join(PROC, "harmonise_log.tsv"), sep="\t", index=False)
    print(log.groupby("outcome")[["n_inst", "n_found", "n_proxy", "n_drop_allele", "n_drop_pal"]].sum())
    print("rows", len(h), "| constants: c_BMI per SD =", round(NH["c_bmi_per_kgm2"] * SD_BMI, 4), "c_ever =", round(C_EVER, 4))
