"""Observational association of ln(PSA) with BMI and smoking in US men (NHANES 2001-2010).

Why: the PSA GWAS (Hoffmann 2025, GCST90461907) adjusted ln(PSA) for BMI and smoking status in MVP and most
other discovery cohorts. A SNP's association with covariate-adjusted PSA is approximately
    b_G,PSA|cov = b_G,PSA - c_BMI * b_G,BMI - c_SMK * b_G,SMK
where c_* are the covariate coefficients of the GWAS regression (Aschard 2015). To reconstruct total effects
we need c_*. MVP's coefficients are not published, so they are estimated here in a comparable population:
men aged >=40 without prostate cancer, PSA 0.01-10 ng/ml (the PSA GWAS inclusion range), model
ln(PSA) ~ age + age^2 + BMI + smoking + race/ethnicity, survey-weighted.

Output: results/tables/nhanes_psa_covariates.tsv and work/nhanes_c.json
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
import statsmodels.api as sm

NH = os.path.join(RAW, "nhanes")
CYC = [("B", "L11PSA_B"), ("C", "L11PSA_C"), ("D", "PSA_D"), ("E", "PSA_E"), ("F", "PSA_F")]


def rd(name):
    return pd.read_sas(os.path.join(NH, name + ".xpt"), format="xport")


frames = []
for s, psaf in CYC:
    psa = rd(psaf)[["SEQN", "LBXP1"]]
    demo = rd("DEMO_" + s)
    wcol = "WTMEC4YR" if s == "B" and "WTMEC4YR" in demo.columns else "WTMEC2YR"
    demo = demo[["SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "SDMVPSU", "SDMVSTRA", "WTMEC2YR"]]
    bmx = rd("BMX_" + s)[["SEQN", "BMXBMI"]]
    smq = rd("SMQ_" + s)[["SEQN", "SMQ020", "SMQ040"]]
    mcq = rd("MCQ_" + s)
    ccols = [c for c in mcq.columns if c.startswith("MCQ230")]
    mcq["pca_hist"] = (mcq[ccols] == 30).any(axis=1)
    mcq = mcq[["SEQN", "MCQ220", "pca_hist"]]
    d = psa.merge(demo, on="SEQN").merge(bmx, on="SEQN", how="left").merge(smq, on="SEQN", how="left").merge(mcq, on="SEQN", how="left")
    d["cycle"] = s
    frames.append(d)
d = pd.concat(frames, ignore_index=True)
n0 = len(d)
flow = [("PSA records 2001-2010", n0)]
d = d[(d.RIAGENDR == 1) & (d.RIDAGEYR >= 40)]
flow.append(("men aged >=40", len(d)))
d = d[d.LBXP1.notna()]
flow.append(("PSA measured", len(d)))
d = d[~d.pca_hist.fillna(False).astype(bool)]
flow.append(("no prostate cancer history", len(d)))
d = d[(d.LBXP1 >= 0.01) & (d.LBXP1 <= 10)]
flow.append(("PSA 0.01-10 ng/ml", len(d)))
d["current_smk"] = np.where(d.SMQ020 == 2, 0, np.where(d.SMQ040.isin([1, 2]), 1, np.where(d.SMQ040 == 3, 0, np.nan)))
d["ever_smk"] = np.where(d.SMQ020 == 1, 1, np.where(d.SMQ020 == 2, 0, np.nan))
d = d.dropna(subset=["BMXBMI", "current_smk", "ever_smk", "WTMEC2YR"])
d = d[d.WTMEC2YR > 0]
flow.append(("complete BMI and smoking", len(d)))
d["lnpsa"] = np.log(d.LBXP1)
d["age"] = d.RIDAGEYR
d["age2"] = (d.age - 60) ** 2
d["w"] = d.WTMEC2YR / len(CYC)
d["psu"] = d.SDMVSTRA.astype(int).astype(str) + "_" + d.SDMVPSU.astype(int).astype(str)
race = pd.get_dummies(d.RIDRETH1.astype(int), prefix="race", drop_first=True).astype(float)


def fit(cols, label):
    X = sm.add_constant(pd.concat([d[cols], race], axis=1).astype(float))
    m = sm.WLS(d.lnpsa, X, weights=d.w).fit(cov_type="cluster", cov_kwds={"groups": d.psu})
    out = []
    for c in cols:
        out.append(dict(model=label, term=c, beta=m.params[c], se=m.bse[c], p=m.pvalues[c], n=int(m.nobs)))
    return out


rows = []
rows += fit(["age", "age2", "BMXBMI", "current_smk"], "BMI + current smoking (primary)")
rows += fit(["age", "age2", "BMXBMI", "ever_smk"], "BMI + ever smoking")
rows += fit(["age", "age2", "BMXBMI"], "BMI only")
rows += fit(["age", "age2", "current_smk"], "current smoking only")
res = pd.DataFrame(rows)
res.to_csv(os.path.join(TAB, "nhanes_psa_covariates.tsv"), sep="\t", index=False)
pd.DataFrame(flow, columns=["step", "n"]).to_csv(os.path.join(TAB, "nhanes_flow.tsv"), sep="\t", index=False)
for f in flow:
    print(*f, sep="\t")
print(res.to_string(index=False))
prim = res[res.model.str.startswith("BMI + current")].set_index("term")
wmean = np.average(d.BMXBMI, weights=d.w)
wsd = math.sqrt(np.average((d.BMXBMI - wmean) ** 2, weights=d.w))
save_json(dict(c_bmi_per_kgm2=prim.loc["BMXBMI", "beta"], c_bmi_se=prim.loc["BMXBMI", "se"],
               c_current_smk=prim.loc["current_smk", "beta"], c_current_smk_se=prim.loc["current_smk", "se"],
               n=int(prim.loc["BMXBMI", "n"]), bmi_mean=wmean, bmi_sd=wsd,
               current_smk_prev=float(np.average(d.current_smk, weights=d.w)),
               median_psa=float(d.LBXP1.median())),
          os.path.join(WORK, "nhanes_c.json"))
print("weighted BMI mean/SD", round(wmean, 2), round(wsd, 2))
