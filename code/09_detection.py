"""Detection-bias decomposition and re-grading of the published prostate-cancer MR claims.

For claim k and prostate-cancer outcome o (primary: Wang 2023 European meta-analysis):
  b_PCa      re-analysed MR estimate (log OR per exposure unit; LD-correlated IVW for drug-target claims, IVW otherwise)
  b_PSA      MR estimate on ln PSA with the BMI/smoking adjustment of the PSA GWAS undone (outcome PSA_total)
  theta_o    MR estimate of ln PSA on prostate cancer in outcome o; primary = cis-KLK instruments (PSA_KLK),
             sensitivity = genome-wide PSA instruments (PSA_gw)
  E          expected log OR produced through PSA alone = b_PSA * theta_o
  b_corr     b_PCa - E (detection-corrected estimate), SE by the delta method (estimates treated as independent)
  DF         detection fraction E / b_PCa, 95% interval from 20,000 Monte Carlo draws (reported when b_PCa replicates)
Re-grade (primary outcome):
  not replicated               b_PCa P >= 0.05, or direction opposite to the published estimate
  not detection-related        replicated, and no PSA pathway in the direction of the association (b_PSA P >= 0.05, or
                               E opposite in sign to b_PCa)
  partly detection-related     replicated, b_PSA P < 0.05, E in the direction of b_PCa, DF < 0.5
  largely detection-related    replicated, b_PSA P < 0.05, E in the direction of b_PCa, DF >= 0.5
DF_adj uses the PSA estimate without undoing the BMI/smoking adjustment (sensitivity).
Outputs: results/tables/detection_<theta>.tsv, results/tables/theta.tsv, results/tables/regrade.tsv
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
from scipy.stats import norm

mr = pd.read_csv(os.path.join(TAB, "mr_all.tsv"), sep="\t")
cm = pd.read_csv(os.path.join(PROC, "claims_map.tsv"), sep="\t")
PCA = ["PCa_Wang2023", "PCa_PRACTICAL", "PCa_UKB", "PCa_FinnGen"]
REGION = {"claim50", "claim51", "claim52", "PSA_KLK"}
rng = np.random.default_rng(20261113)


def primary(key, outcome):
    s = mr[(mr.key == key) & (mr.outcome == outcome)]
    if s.empty:
        return None
    for m in (["IVW (LD-correlated)"] if key in REGION else []) + ["IVW", "Wald ratio"]:
        t = s[s.method == m]
        if len(t):
            return t.iloc[0]
    return None


# theta per outcome
th_rows = []
for key in ("PSA_KLK", "PSA_gw"):
    for o in PCA + ["PSAtest_UKB"]:
        s = mr[(mr.key == key) & (mr.outcome == o)]
        for r in s.itertuples():
            th_rows.append(dict(instruments=key, outcome=o, method=r.method, nsnp=r.nsnp, theta=r.b, se=r.se, p=r.p,
                                or_per_lnpsa=math.exp(r.b), lo=math.exp(r.b - 1.96 * r.se), hi=math.exp(r.b + 1.96 * r.se)))
theta = pd.DataFrame(th_rows)
theta.to_csv(os.path.join(TAB, "theta.tsv"), sep="\t", index=False)
print(theta.to_string(index=False))


def decompose(thkey):
    rows = []
    for c in cm[cm.inst_mode != "exclude"].itertuples():
        key = f"claim{c.claim:02d}"
        psa = primary(key, "PSA_total")
        psa_adj = primary(key, "PSA_adj")
        test = primary(key, "PSAtest_UKB")
        for o in PCA:
            y = primary(key, o)
            t = primary(thkey, o)
            if y is None or t is None or psa is None:
                continue
            E = psa.b * t.b
            seE = math.sqrt((psa.se * t.b) ** 2 + (psa.b * t.se) ** 2)
            bc = y.b - E
            sc = math.sqrt(y.se ** 2 + seE ** 2)
            dr = rng.normal(y.b, y.se, 20000)
            de = rng.normal(psa.b, psa.se, 20000) * rng.normal(t.b, t.se, 20000)
            df = de / dr
            rows.append(dict(claim=c.claim, exposure=c.exposure, category=c.category, grade=c.grade, orig_or=c.orig_or,
                             orig_outcome_type=c.orig_outcome_type, outcome=o, method=y.method, nsnp=int(y.nsnp),
                             b_pca=y.b, se_pca=y.se, p_pca=y.p, b_psa=psa.b, se_psa=psa.se, p_psa=psa.p,
                             b_psa_adj=psa_adj.b if psa_adj is not None else np.nan,
                             p_psa_adj=psa_adj.p if psa_adj is not None else np.nan,
                             b_test=test.b if test is not None else np.nan, se_test=test.se if test is not None else np.nan,
                             p_test=test.p if test is not None else np.nan,
                             theta=t.b, se_theta=t.se, E=E, se_E=seE, b_corr=bc, se_corr=sc, p_corr=2 * norm.sf(abs(bc / sc)),
                             DF=E / y.b, DF_lo=np.percentile(df, 2.5), DF_hi=np.percentile(df, 97.5),
                             DF_adj=(psa_adj.b * t.b / y.b) if psa_adj is not None else np.nan))
    return pd.DataFrame(rows)


def regrade(d):
    d = d[d.outcome == "PCa_Wang2023"].copy()
    same_dir = np.sign(d.b_pca) == np.sign(np.log(d.orig_or))
    rep = (d.p_pca < 0.05) & same_dir
    det = (d.p_psa < 0.05) & (np.sign(d.E) == np.sign(d.b_pca))
    d["replicated"] = rep
    d["new_grade"] = np.where(~rep, "not replicated",
                              np.where(~det, "not detection-related",
                                       np.where(d.DF >= 0.5, "largely detection-related", "partly detection-related")))
    return d


for thkey in ("PSA_KLK", "PSA_gw"):
    d = decompose(thkey)
    d.to_csv(os.path.join(TAB, f"detection_{thkey}.tsv"), sep="\t", index=False)
    g = regrade(d)
    g.to_csv(os.path.join(TAB, f"regrade_{thkey}.tsv"), sep="\t", index=False)
    print(f"\n== theta from {thkey} ==")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(pd.crosstab(g.grade, g.new_grade, margins=True))
    show = g[["claim", "exposure", "grade", "orig_or", "b_pca", "p_pca", "b_psa", "p_psa", "E", "DF", "DF_lo", "DF_hi", "DF_adj",
              "p_corr", "p_test", "new_grade"]]
    print(show.round(4).to_string(index=False))
