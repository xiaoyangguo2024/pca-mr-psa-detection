"""How much prostate cancer 'effect' can a PSA-only exposure produce? Analytic model + MR simulation.

Model (one follow-up period, men 55-69):
  benign ln PSA ~ N(mu0, sd0)  (NHANES 2001-2010, men 55-69 without prostate cancer, PSA 0.01-10 not truncated here)
  latent cancer with prevalence pi; cancer adds delta_c to ln PSA
  a man is PSA-tested with probability s (screening intensity); PSA > T leads to biopsy with probability u;
  biopsy finds a latent cancer with sensitivity b; latent cancers are also found clinically with probability k,
  independent of PSA.  Detected cancer = case; everyone else = control.
  An exposure shifts ln PSA by dpsa per unit and has NO effect on cancer.
Analytic part: theta_sim(s, T) = d logOR(detected PCa) / d ln PSA, evaluated at dpsa = 0.
MR part: 20 SNPs -> exposure X -> ln PSA only; two-sample MR (IVW) of X on detected PCa, 200 replicates per scenario.
Outputs: results/tables/sim_theta.tsv, results/tables/sim_mr.tsv
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
from scipy.stats import norm

NH = os.path.join(RAW, "nhanes")
CYC = [("B", "L11PSA_B"), ("C", "L11PSA_C"), ("D", "PSA_D"), ("E", "PSA_E"), ("F", "PSA_F")]
fr = []
for s, pf in CYC:
    psa = pd.read_sas(os.path.join(NH, pf + ".xpt"), format="xport")[["SEQN", "LBXP1"]]
    demo = pd.read_sas(os.path.join(NH, "DEMO_" + s + ".xpt"), format="xport")[["SEQN", "RIAGENDR", "RIDAGEYR", "WTMEC2YR"]]
    mcq = pd.read_sas(os.path.join(NH, "MCQ_" + s + ".xpt"), format="xport")
    cc = [c for c in mcq.columns if c.startswith("MCQ230")]
    mcq["pca"] = (mcq[cc] == 30).any(axis=1)
    fr.append(psa.merge(demo, on="SEQN").merge(mcq[["SEQN", "pca"]], on="SEQN", how="left"))
d = pd.concat(fr)
d = d[(d.RIAGENDR == 1) & d.RIDAGEYR.between(55, 69) & d.LBXP1.notna() & (d.LBXP1 > 0) & ~d.pca.fillna(False).astype(bool)
      & (d.WTMEC2YR > 0)]
lp = np.log(d.LBXP1)
MU0 = float(np.average(lp, weights=d.WTMEC2YR))
SD0 = float(math.sqrt(np.average((lp - MU0) ** 2, weights=d.WTMEC2YR)))
print(f"NHANES men 55-69 n={len(d)}: ln PSA mean {MU0:.3f} sd {SD0:.3f} (median PSA {d.LBXP1.median():.2f})")

BASE = dict(pi=0.25, delta_c=0.5, u=0.5, b=0.7, k=0.02)


def p_detect(shift, s, T, pi, delta_c, u, b, k):
    p_hi = 1 - norm.cdf((math.log(T) - MU0 - delta_c - shift) / SD0)
    det = 1 - (1 - s * u * b * p_hi) * (1 - k)
    return pi * det


def theta(s, T, h=1e-3, **kw):
    par = {**BASE, **kw}
    r1, r0 = p_detect(h, s, T, **par), p_detect(-h, s, T, **par)
    return (math.log(r1 / (1 - r1)) - math.log(r0 / (1 - r0))) / (2 * h)


rows = []
for T in (3.0, 4.0):
    for s in np.round(np.arange(0, 1.0001, 0.05), 2):
        rows.append(dict(scenario="base", T=T, s=s, theta=theta(s, T), risk=p_detect(0, s, T, **BASE)))
    # sensitivity on the latent-cancer parameters (s = 0.3 and 0.7)
    for s in (0.3, 0.7):
        for name, vals in (("pi", (0.15, 0.35)), ("delta_c", (0.25, 0.75)), ("u", (0.3, 0.8)), ("b", (0.5, 0.9)), ("k", (0.01, 0.05))):
            for v in vals:
                rows.append(dict(scenario=f"{name}={v}", T=T, s=s, theta=theta(s, T, **{name: v}),
                                 risk=p_detect(0, s, T, **{**BASE, name: v})))
th = pd.DataFrame(rows)
th.to_csv(os.path.join(TAB, "sim_theta.tsv"), sep="\t", index=False)
print(th[th.scenario == "base"].pivot(index="s", columns="T", values="theta").iloc[::4].round(3))

# ---- MR simulation -------------------------------------------------------------------------------------------------
rng = np.random.default_rng(20261113)
NSNP, N = 20, 200_000
maf = rng.uniform(0.1, 0.5, NSNP)
gx = rng.normal(0, 1, NSNP); gx = gx / np.sqrt(np.sum(gx ** 2 * 2 * maf * (1 - maf))) * math.sqrt(0.05)  # R2_X = 5%


def one(dpsa, s, T, rep):
    r = np.random.default_rng(rep)
    # sample 1: exposure GWAS
    G1 = r.binomial(2, maf, size=(N, NSNP)).astype(np.float32)
    X1 = G1 @ gx + r.normal(0, math.sqrt(0.95), N)
    Gc = G1 - G1.mean(0)
    bx = (Gc * (X1 - X1.mean())[:, None]).sum(0) / (Gc ** 2).sum(0)
    sx = np.sqrt(np.var(X1) / (Gc ** 2).sum(0))
    # sample 2: outcome GWAS (detected prostate cancer)
    G2 = r.binomial(2, maf, size=(N, NSNP)).astype(np.float32)
    X2 = G2 @ gx + r.normal(0, math.sqrt(0.95), N)
    can = r.random(N) < BASE["pi"]
    lpsa = MU0 + dpsa * X2 + BASE["delta_c"] * can + r.normal(0, SD0, N)
    tested = r.random(N) < s
    biop = tested & (lpsa > math.log(T)) & (r.random(N) < BASE["u"])
    det = can & ((biop & (r.random(N) < BASE["b"])) | (r.random(N) < BASE["k"]))
    y = det
    # per-allele log OR from allele counts
    ca, co = G2[y].sum(0), G2[~y].sum(0)
    na, no = 2 * y.sum(), 2 * (~y).sum()
    by = np.log(ca / (na - ca)) - np.log(co / (no - co))
    sy = np.sqrt(1 / ca + 1 / (na - ca) + 1 / co + 1 / (no - co))
    w = bx ** 2 / sy ** 2
    ivw = np.sum(bx * by / sy ** 2) / np.sum(bx ** 2 / sy ** 2)
    se = 1 / math.sqrt(np.sum(w))
    return ivw, se, y.mean()


mrrows = []
for dpsa in (-0.10, -0.05):
    for s in (0.1, 0.3, 0.5, 0.7, 0.9):
        T = 4.0
        res = np.array([one(dpsa, s, T, 1000 * int(s * 10) + rep + (0 if dpsa < -0.07 else 500000)) for rep in range(200)])
        exp_th = theta(s, T) * dpsa
        mrrows.append(dict(dpsa=dpsa, s=s, T=T, ivw_mean=res[:, 0].mean(), ivw_sd=res[:, 0].std(), se_mean=res[:, 1].mean(),
                           power=float(np.mean(np.abs(res[:, 0] / res[:, 1]) > 1.96)), expected=exp_th,
                           case_frac=res[:, 2].mean()))
        print(mrrows[-1], flush=True)
pd.DataFrame(mrrows).to_csv(os.path.join(TAB, "sim_mr.tsv"), sep="\t", index=False)
