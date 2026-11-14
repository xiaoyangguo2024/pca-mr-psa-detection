"""Manuscript tables and the supplementary workbook; also prints the numbers quoted in the text.

Outputs
  results/tables/Table1_sources.tsv        data sources
  results/tables/Table2_replicated.tsv     detection decomposition for replicated claims (cis-KLK theta)
  results/Supplementary_Tables.xlsx        S1-S10
  results/tables/summary_numbers.json      numbers quoted in the abstract and results
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
from scipy.stats import spearmanr

cm = pd.read_csv(os.path.join(PROC, "claims_map.tsv"), sep="\t")
g = pd.read_csv(os.path.join(TAB, "regrade_PSA_KLK.tsv"), sep="\t")
g2 = pd.read_csv(os.path.join(TAB, "regrade_PSA_gw.tsv"), sep="\t")
d = pd.read_csv(os.path.join(TAB, "detection_PSA_KLK.tsv"), sep="\t")
dg = pd.read_csv(os.path.join(TAB, "detection_PSA_gw.tsv"), sep="\t")
mr = pd.read_csv(os.path.join(TAB, "mr_all.tsv"), sep="\t")
th = pd.read_csv(os.path.join(TAB, "theta.tsv"), sep="\t")
isum = pd.read_csv(os.path.join(PROC, "instrument_summary.tsv"), sep="\t")
psa = pd.read_csv(os.path.join(PROC, "psa_instruments.tsv"), sep="\t")
hlog = pd.read_csv(os.path.join(PROC, "harmonise_log.tsv"), sep="\t")

# exposure scale of each exposure GWAS (effects are per one unit of this scale)
UNIT = {1: "log odds (linear-scale GWAS)", 2: "SD", 3: "log odds", 4: "SD (years of schooling)", 5: "SD", 6: "years",
        7: "SD", 8: "SD", 10: "SD", 11: "SD", 12: "SD", 13: "SD", 14: "SD", 15: "SD", 16: "SD", 17: "SD", 19: "SD",
        20: "SD", 21: "SD", 22: "SD", 23: "SD", 24: "SD", 25: "SD", 27: "SD", 28: "SD", 29: "SD", 30: "SD", 31: "SD",
        32: "SD", 33: "SD", 34: "SD", 35: "SD", 36: "SD", 37: "SD", 38: "SD", 39: "SD", 40: "SD", 41: "SD", 42: "SD",
        43: "SD", 45: "SD", 46: "SD", 47: "SD", 48: "SD", 50: "SD LDL-C", 51: "mmHg SBP", 52: "SD HbA1c", 53: "log odds",
        54: "log odds", 55: "log odds"}


def ci(b, se, f=np.exp, nd=2):
    return f"{f(b):.{nd}f} ({f(b - 1.96 * se):.{nd}f}-{f(b + 1.96 * se):.{nd}f})"


def pfmt(p):
    return f"{p:.1e}" if p < 0.001 else f"{p:.3f}"


# ---------------------------------------------------------------- Table 2
r = g[g.replicated].sort_values("DF", ascending=False)
t2 = pd.DataFrame({
    "Exposure (unit)": [f"{nice(x.exposure)} ({UNIT.get(x.claim, 'unit')})" for x in r.itertuples()],
    "Original grade": r.grade,
    "Original OR": r.orig_or,
    "Reanalyzed OR (95% CI)": [ci(x.b_pca, x.se_pca) for x in r.itertuples()],
    "P": [pfmt(p) for p in r.p_pca],
    "PSA ratio per unit (95% CI)": [ci(x.b_psa, x.se_psa, nd=3) for x in r.itertuples()],
    "P (PSA)": [pfmt(p) for p in r.p_psa],
    "PSA-driven OR": [f"{math.exp(x):.3f}" for x in r.E],
    "Detection fraction (95% interval)": [f"{x.DF:.2f} ({x.DF_lo:.2f} to {x.DF_hi:.2f})" for x in r.itertuples()],
    "Corrected OR (95% CI)": [ci(x.b_corr, x.se_corr) for x in r.itertuples()],
    "P (corrected)": [pfmt(p) for p in r.p_corr],
    "Regrade": r.new_grade,
})
t2.to_csv(os.path.join(TAB, "Table2_replicated.tsv"), sep="\t", index=False)

# ---------------------------------------------------------------- Table 1
t1 = pd.DataFrame([
    ("Prostate cancer (primary)", "Wang et al. 2023, European meta-analysis (PRACTICAL/ELLIPSE, UK Biobank, FinnGen, MVP, US biobanks)", "122,188 / 604,640", "GCST90274714", "log OR"),
    ("Prostate cancer", "Schumacher et al. 2018, PRACTICAL/OncoArray", "79,148 / 61,106", "ieu-b-85", "log OR"),
    ("Prostate cancer", "UK Biobank (MRC-IEU OpenGWAS)", "9,132 / 173,493", "ieu-b-4809", "linear model; converted to log OR"),
    ("Prostate cancer", "FinnGen R5, C3_PROSTATE_EXALLC", "6,311 / 74,685", "finn-b-C3_PROSTATE_EXALLC", "log OR"),
    ("ln PSA", "Hoffmann et al. 2025, multi-ancestry (76% European; 96.5% MVP)", "392,522 men without prostate cancer", "GCST90461907", "ln(ng/ml); adjusted for BMI and smoking in most cohorts"),
    ("Ever had a PSA test", "UK Biobank field 2365 (MRC-IEU pipeline)", "62,489 / 137,921", "ukb-b-15695", "linear model; converted to log OR"),
    ("BMI (correction of PSA effects)", "UK Biobank (MRC-IEU pipeline)", "461,460", "ukb-b-19953", "SD (4.8 kg/m²)"),
    ("Smoking initiation (correction of PSA effects)", "GSCAN (Liu et al. 2019), excluding 23andMe", "632,802", "ieu-b-4877", "log odds"),
    ("PSA covariate coefficients", "NHANES 2001-2010, men ≥40 y without prostate cancer, PSA 0.01-10 ng/ml", f"{load_json(os.path.join(WORK, 'nhanes_c.json'))['n']:,}", "CDC NHANES", "survey-weighted regression"),
], columns=["Use", "Source", "Cases / controls or N", "Accession", "Scale"])
t1.to_csv(os.path.join(TAB, "Table1_sources.tsv"), sep="\t", index=False)

# ---------------------------------------------------------------- supplementary workbook
s1 = cm.merge(isum.rename(columns={"key": "inst_key"}).assign(claim=lambda x: pd.to_numeric(x.inst_key.str.replace("claim", ""), errors="coerce")),
              on="claim", how="left", suffixes=("", "_inst"))
s1 = s1[["claim", "category", "exposure", "source", "source_pmid", "grade", "orig_or", "orig_lo", "orig_hi", "orig_unit",
         "orig_outcome", "orig_outcome_type", "orig_n_iv", "exp_id", "exp_trait", "exp_author", "exp_n", "inst_mode",
         "p_used", "n_snp", "F_mean", "F_min", "region_genes", "map_note"]]
s1["exp_unit"] = s1.claim.map(UNIT)
sheets = {
    "S1_claims_instruments": s1,
    "S2_PSA_instruments": psa,
    "S3_theta": th,
    "S4_MR_all": mr,
    "S5_detection_KLK": d,
    "S6_detection_genomewide": dg,
    "S7_regrade": g.merge(g2[["claim", "DF", "new_grade"]].rename(columns={"DF": "DF_gw", "new_grade": "new_grade_gw"}), on="claim"),
    "S8_NHANES": pd.read_csv(os.path.join(TAB, "nhanes_psa_covariates.tsv"), sep="\t"),
    "S9_simulation_theta": pd.read_csv(os.path.join(TAB, "sim_theta.tsv"), sep="\t"),
    "S10_simulation_MR": pd.read_csv(os.path.join(TAB, "sim_mr.tsv"), sep="\t"),
    "S11_harmonisation": hlog,
}
# detection fraction with each estimator: exposure effects on PSA and on prostate cancer from the same estimator,
# theta fixed at the primary cis-KLK LD-aware IVW estimate in the Wang 2023 meta-analysis
th0 = float(th[(th.instruments == "PSA_KLK") & (th.outcome == "PCa_Wang2023") & (th.method == "IVW (LD-correlated)")].theta.iloc[0])
est = []
for x in g[g.replicated].itertuples():
    k = f"claim{x.claim:02d}"
    for m in ("IVW", "Weighted median", "Weighted mode", "MR-Egger"):
        a = mr[(mr.key == k) & (mr.outcome == "PCa_Wang2023") & (mr.method == m)]
        b = mr[(mr.key == k) & (mr.outcome == "PSA_total") & (mr.method == m)]
        if len(a) and len(b):
            est.append(dict(claim=x.claim, exposure=x.exposure, re_grade=x.new_grade, method=m, nsnp=int(a.nsnp.iloc[0]),
                            b_pca=a.b.iloc[0], se_pca=a.se.iloc[0], p_pca=a.p.iloc[0], b_psa=b.b.iloc[0],
                            se_psa=b.se.iloc[0], p_psa=b.p.iloc[0], egger_int_p_pca=a.egger_int_p.iloc[0],
                            egger_int_p_psa=b.egger_int_p.iloc[0], DF=b.b.iloc[0] * th0 / a.b.iloc[0]))
dfm = pd.DataFrame(est)
dfm.to_csv(os.path.join(TAB, "DF_by_estimator.tsv"), sep="\t", index=False)
sheets["S12_DF_by_estimator"] = dfm
mv = os.path.join(TAB, "mvmr.tsv")
if os.path.exists(mv):
    sheets["S13_MVMR"] = pd.read_csv(mv, sep="\t")
with pd.ExcelWriter(os.path.join(RES, "Supplementary_Tables.xlsx")) as xw:
    for k, v in sheets.items():
        v.to_excel(xw, sheet_name=k[:31], index=False)

# ---------------------------------------------------------------- numbers for the text
same = cm.map_note.str.startswith("Same")
n = {}
n["n_claims"], n["n_reanalysed"] = len(cm), int((cm.inst_mode != "exclude").sum())
n["n_replicated"] = int(g.replicated.sum())
n["replicated_by_grade"] = g.groupby("grade").replicated.agg(["sum", "size"]).astype(int).to_dict("index")
gg = g.merge(cm[["claim", "map_note"]], on="claim")
n["replicated_same_gwas"] = [int(gg[gg.map_note.str.startswith("Same")].replicated.sum()), int(gg.map_note.str.startswith("Same").sum())]
n["replicated_overall_outcome"] = [int(gg[gg.orig_outcome_type == "overall"].replicated.sum()), int((gg.orig_outcome_type == "overall").sum())]
n["new_grade_KLK"] = g.new_grade.value_counts().to_dict()
n["new_grade_gw"] = g2.new_grade.value_counts().to_dict()
det = g[g.new_grade.isin(["largely detection-related", "partly detection-related"])]
n["DF_median_KLK"], n["DF_range_KLK"] = float(det.DF.median()), [float(det.DF.min()), float(det.DF.max())]
rr = g[(g.grade == "Robust") & g.replicated]
n["robust_replicated_detection"] = [int(rr.claim.isin(det.claim).sum()), len(rr)]
n["corrected_nonsignificant"] = list(det[det.p_corr >= 0.05].exposure)
z_obs, z_pred = g.b_pca / g.se_pca, g.E / g.se_pca
rho = spearmanr(z_obs, z_pred)
n["spearman_obs_pred_all"] = [float(rho.statistic), float(rho.pvalue)]
nr = g[~g.replicated]
rho2 = spearmanr(nr.b_pca / nr.se_pca, nr.E / nr.se_pca)
n["spearman_obs_pred_notreplicated"] = [float(rho2.statistic), float(rho2.pvalue)]
k = th[(th.instruments == "PSA_KLK") & (th.method == "IVW (LD-correlated)")].set_index("outcome")
n["theta_KLK"] = {o: [round(float(k.loc[o, "theta"]), 3), round(float(k.loc[o, "se"]), 3)] for o in k.index}
n["n_psa_instruments"] = psa.key.value_counts().to_dict()
n["bmi"] = g[g.claim == 8][["b_pca", "se_pca", "b_psa", "se_psa", "b_psa_adj", "DF", "DF_lo", "DF_hi", "DF_adj", "b_corr", "se_corr", "p_corr"]].iloc[0].to_dict()
if os.path.exists(mv):
    m = pd.read_csv(mv, sep="\t")
    m = m[(m.method == "MV-IVW") & (m.term == "X")].assign(claim=lambda x: x.key.str.replace("claim", "").astype(int))
    m = m.merge(g[["claim", "exposure", "b_pca", "se_pca", "new_grade"]], on="claim")
    m["ratio"], m["se_ratio"] = m.b / m.b_pca, m.se / m.se_pca
    m["detection"] = m.new_grade.isin(["largely detection-related", "partly detection-related"])
    n["mvmr_ratio_median"] = {("detection" if k else "no_psa_pathway"): float(v) for k, v in m.groupby("detection").ratio.median().items()}
    n["mvmr_se_ratio_range"] = [float(m.se_ratio.min()), float(m.se_ratio.max())]
    n["mvmr_n_sig"] = int((m.p < 0.05).sum())
    n["mvmr_lnpsa_range"] = [float(x) for x in pd.read_csv(mv, sep="\t").query("method == 'MV-IVW' and term == 'lnPSA'").b.agg(["min", "max"])]
    print(m[["exposure", "b", "se", "p", "b_pca", "ratio", "se_ratio", "new_grade"]].round(3).to_string(index=False))
save_json(n, os.path.join(TAB, "summary_numbers.json"))
print(json.dumps(n, indent=1, default=float))
print(t2.to_string(index=False))
