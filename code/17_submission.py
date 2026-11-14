"""Assemble the submission package for The Prostate (Wiley; Research Exchange) in submission/.

Contents (named as the journal asks: one file per figure and per table, named by number)
  Manuscript.docx                                   title page with required statements, abstract, text, references,
                                                    figure legends, Supporting Information list (16_build_docx.py)
  Table 1.docx, Table 2.docx                        Word tables, one per file
  Figure 1.tif ... Figure 4.tif                     TIFF, 1200 dpi, RGB, LZW (13_figures.py)
  Supporting Information - Figures S1-S2.pdf
  Supporting Information - Tables S1-S13.xlsx
  Supporting Information - STROBE-MR checklist.docx  items paraphrased, with manuscript locations
  Cover letter.docx                                 author fields left as placeholders
Only files to upload go into submission/. Elsewhere:
  code_release/code_and_derived_tables.zip          code/, derived tables and README for the public repository; never
                                                    the OpenGWAS token, env/, Rlib/, data/raw/ or data/ref/
  manuscript/manuscript_with_figures.docx           reading copy for the authors (16_build_docx.py)
Run: python3 -I code/17_submission.py
"""
import sys, os, shutil, zipfile; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
import fitz
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_LINE_SPACING
from docx.oxml.ns import qn

SUB = os.path.join(ROOT, "submission")
MS = os.path.join(ROOT, "manuscript")
if os.path.isdir(SUB):
    shutil.rmtree(SUB)
os.makedirs(SUB)
NUM = load_json(os.path.join(TAB, "summary_numbers.json"))
TITLE = open(os.path.join(MS, "manuscript.md"), encoding="utf-8").readline().lstrip("# ").strip()


def new_doc(size=11):
    d = Document(); st = d.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(size)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    st.paragraph_format.space_after = Pt(6); st.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
    for s in d.sections:
        s.page_width, s.page_height = Cm(21.0), Cm(29.7)
        s.left_margin = s.right_margin = s.top_margin = s.bottom_margin = Cm(2.5)
    return d


# ---------------------------------------------------------------- manuscript, tables, figures
shutil.copy(os.path.join(MS, "manuscript.docx"), os.path.join(SUB, "Manuscript.docx"))
for k in (1, 2):
    shutil.copy(os.path.join(MS, f"Table {k}.docx"), os.path.join(SUB, f"Table {k}.docx"))
for k, f in enumerate(["Fig1_concept_flow", "Fig2_theta", "Fig3_detection_fraction", "Fig4_regrade"], start=1):
    shutil.copy(os.path.join(FIG, f + ".tif"), os.path.join(SUB, f"Figure {k}.tif"))

# ---------------------------------------------------------------- cover letter
d = new_doc(11)
for t in ("[Date]", "", "The Editor", "The Prostate", "", "Dear Editor,"):
    d.add_paragraph(t)
d.add_paragraph(f"We submit our manuscript entitled “{TITLE}” for consideration as an Original Article in "
                "The Prostate.")
d.add_paragraph(
    "Mendelian randomization (MR) studies are now the most common source of causal claims about risk factors for "
    "prostate cancer. Every such study uses diagnosed prostate cancer as its outcome, and diagnosis depends heavily on "
    "PSA testing. An exposure that changes serum PSA therefore changes the chance that a latent cancer is found, even if "
    "it has no effect on cancer development. To our knowledge, no previous study has quantified this detection pathway "
    "across the published MR evidence.")
d.add_paragraph(
    "We reanalyzed 50 of the 55 statistically significant MR claims graded in a 2024 umbrella review, using public "
    "summary statistics, including a GWAS of PSA levels in 392,522 men and a prostate cancer meta-analysis of 122,188 "
    "cases. Using KLK3 variants, we estimated that each unit of genetically higher ln PSA doubles the odds of a prostate "
    f"cancer diagnosis. Of the {NUM['n_replicated']} claims that replicated, 10 were partly explained by the exposure's "
    "effect on PSA, including every replicated claim graded robust. For body mass index, about three-quarters of the "
    "association was attributable to detection, and the corrected estimate was null. Associations for smoking, "
    "education, and testosterone were also no longer statistically significant after correction, whereas the "
    "association with telomere length persisted.")
d.add_paragraph(
    "We believe these findings will interest the readers of The Prostate who design, review, or cite MR studies of "
    "prostate cancer. We propose a simple check, using public data, that can be added to every such study. All "
    "analyses use publicly available data, and the code is provided.")
d.add_paragraph(
    "This manuscript has not been published and is not under consideration elsewhere. All authors have approved the "
    "submission and declare [no conflicts of interest]. [Statement on the use of AI tools, to be confirmed by the "
    "authors.]")
d.add_paragraph("Thank you for considering our work.")
d.add_paragraph("Sincerely,")
d.add_paragraph("[Corresponding author name, affiliation, postal address, e-mail]")
d.save(os.path.join(SUB, "Cover letter.docx"))

# ---------------------------------------------------------------- Supporting Information: figures
LEG = {
    "FigS1_cohorts": ("Figure S1.", "Estimates for the 16 replicated claims in the Wang et al meta-analysis, PRACTICAL, "
                      "UK Biobank, and FinnGen release 5, each divided by the absolute Wang et al estimate (95% CI; "
                      "intervals clipped at ±4). PRACTICAL, UK Biobank, and FinnGen are components of the "
                      "meta-analysis. ALT, alanine aminotransferase; BMI, body mass index; CCL2, C-C motif chemokine "
                      "ligand 2; CI, confidence interval; IGF-I, insulin-like growth factor I; LDL, low-density "
                      "lipoprotein."),
    "FigS2_observed_vs_predicted": ("Figure S2.", "Observed log odds ratio against the log odds ratio predicted through "
                                    "PSA (b_PSA × θ) for all 50 reanalyzed claims, both divided by the "
                                    "standard error of the observed estimate. The dashed line is the line of identity. "
                                    "PSA, prostate-specific antigen; θ, log odds ratio of diagnosis per unit ln PSA."),
}
out = fitz.open()
p = out.new_page(width=595, height=842)
p.insert_htmlbox(fitz.Rect(60, 80, 535, 400),
                 f"<p style='font-family:serif;font-size:13pt'><b>Supporting Information</b><br><br>{TITLE}<br><br>"
                 "Figures S1–S2</p>")
for f, (lab, leg) in LEG.items():
    src = fitz.open(os.path.join(FIG, f + ".pdf"))
    r = src[0].rect
    page = out.new_page(width=595, height=842)
    w = 475; h = w * r.height / r.width
    if h > 600:
        h = 600; w = h * r.width / r.height
    x0 = (595 - w) / 2
    page.show_pdf_page(fitz.Rect(x0, 60, x0 + w, 60 + h), src, 0)
    page.insert_htmlbox(fitz.Rect(60, 75 + h, 535, 820),
                        f"<p style='font-family:serif;font-size:10pt'><b>{lab}</b> {leg}</p>")
out.save(os.path.join(SUB, "Supporting Information - Figures S1-S2.pdf"))

# ---------------------------------------------------------------- Supporting Information: tables
shutil.copy(os.path.join(RES, "Supplementary_Tables.xlsx"), os.path.join(SUB, "Supporting Information - Tables S1-S13.xlsx"))

# ---------------------------------------------------------------- Supporting Information: STROBE-MR checklist
ITEMS = [
    ("Title and abstract", "1", "MR design indicated in title and/or abstract", "Title; Abstract"),
    ("Background", "2", "Scientific background and rationale", "Introduction, paragraphs 1–4"),
    ("Objectives", "3", "Objectives, including the causal question", "Introduction, last paragraph"),
    ("Study design and data sources", "4a", "Study design and underlying population", "Materials and Methods: Design and claims; Outcomes"),
    ("", "4b", "Data sources: each GWAS, sample size, ancestry, covariates", "Materials and Methods: Outcomes; Table 1; Table S1"),
    ("", "4c", "Ethics approval and consent", "Title page: ethics approval and patient consent statements"),
    ("Assumptions", "5", "Assumptions of the instrumental-variable analysis", "Materials and Methods: Effect of PSA on diagnosis; Statistical analysis; Discussion (limitations)"),
    ("Statistical methods: main analysis", "6a", "Selection of genetic variants", "Materials and Methods: Genetic instruments; Effect of PSA on diagnosis"),
    ("", "6b", "Harmonization and handling of variants", "Materials and Methods: Outcomes (last paragraph); Table S11"),
    ("", "6c", "MR estimator and how estimates were combined", "Materials and Methods: Statistical analysis"),
    ("", "6d", "Scale of estimates; conversions", "Materials and Methods: Outcomes; Undoing the BMI and smoking adjustment"),
    ("", "6e", "Multiple testing", "Materials and Methods: Statistical analysis"),
    ("Assessment of assumptions", "7", "Methods to assess the assumptions", "Materials and Methods: Statistical analysis (MR-Egger, weighted median, weighted mode, Q)"),
    ("Sensitivity and additional analyses", "8", "Sensitivity and additional analyses", "Materials and Methods: Statistical analysis; Analytic model and simulation"),
    ("Software and preregistration", "9", "Software and versions; preregistration", "Materials and Methods: Statistical analysis"),
    ("Descriptive data", "10", "Characteristics of data sources and instruments", "Table 1; Tables S1–S2"),
    ("Main results", "11", "Main MR results with precision", "Results; Table 2; Figures 2–4"),
    ("Assessment of assumptions", "12", "Results of assumption checks", "Results: Sensitivity analyses; Tables S4, S12"),
    ("Sensitivity and additional analyses", "13", "Results of sensitivity and additional analyses", "Results: Sensitivity analyses; Tables S6–S7, S9–S10, S12–S13"),
    ("Key results", "14", "Summary of key results", "Discussion, paragraphs 1–2"),
    ("Limitations", "15", "Limitations, including violations of assumptions", "Discussion (limitations)"),
    ("Interpretation", "16", "Interpretation and implications", "Discussion; Conclusions"),
    ("Generalizability", "17", "Generalizability", "Discussion (limitations: generalizability)"),
    ("Funding", "18", "Sources of funding", "Title page: funding statement"),
    ("Data and data sharing", "19", "Availability of data and code", "Title page: data availability statement"),
    ("Conflicts of interest", "20", "Conflicts of interest", "Title page: conflict of interest disclosure"),
]
d = new_doc(10)
d.add_paragraph().add_run("Supporting Information: STROBE-MR checklist").bold = True
d.add_paragraph("Items are summarized from the STROBE-MR statement (Skrivankova et al, JAMA 2021); see the statement for "
                "the full wording.")
tb = d.add_table(rows=1, cols=4); tb.style = "Table Grid"
for j, h in enumerate(("Section", "Item", "Topic", "Location in manuscript")):
    tb.cell(0, j).text = ""; r = tb.cell(0, j).paragraphs[0].add_run(h); r.bold = True
for it in ITEMS:
    row = tb.add_row().cells
    for j, v in enumerate(it):
        row[j].text = v
d.save(os.path.join(SUB, "Supporting Information - STROBE-MR checklist.docx"))

# ---------------------------------------------------------------- code and derived tables
README = f"""# {TITLE}

Code and derived tables. All inputs are public: GWAS Catalog (GCST90274714, GCST90461907), IEU OpenGWAS and NHANES.

Steps (run from the project root; Python 3.9 with pandas, numpy, scipy, statsmodels, matplotlib, python-docx, PyMuPDF;
R 4.5.3 with MendelianRandomization 0.10.0 and data.table):
  01_download_outcomes.sh / 01b_download_one.sh  GWAS Catalog harmonized files (MD5 checked)
  03_nhanes_psa_covariates.py                     NHANES coefficients of ln PSA on BMI and smoking
  04_claims_map.py                                claims from Cui et al 2024 mapped to exposure GWAS
  05_instruments.py, 06b_psa_instruments.py       instruments (OpenGWAS clumping)
  06_lookups.py                                   instrument lookups in outcome GWAS
  07_harmonise.py                                 harmonization; BMI/smoking adjustment of the PSA GWAS undone
  08_mr.R                                         univariable MR
  09_detection.py                                 theta, detection decomposition and regrading
  10_simulation.py                                analytic model and MR simulation
  11_mvmr_prep.py, 12_mvmr.R                      multivariable MR
  13_figures.py, 14_tables.py                     figures and tables
  15_references.py, 16_build_docx.py              reference list (PubMed) and manuscript
  17_submission.py                                submission package

OpenGWAS API calls need a personal token (https://api.opengwas.io). Scripts read it from a file named
.opengwas_token in the project root; the token is not included here.
"""
REL = os.path.join(ROOT, "code_release")
os.makedirs(REL, exist_ok=True)
zp = os.path.join(REL, "code_and_derived_tables.zip")
with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
    z.writestr("README.md", README)
    for f in sorted(os.listdir(os.path.join(ROOT, "code"))):
        if f.endswith((".py", ".R", ".sh")):
            z.write(os.path.join(ROOT, "code", f), f"code/{f}")
    for f in sorted(os.listdir(TAB)):
        if f.endswith((".tsv", ".json")):
            z.write(os.path.join(TAB, f), f"results/tables/{f}")
    for f in ("claims_map.tsv", "instrument_summary.tsv", "psa_instruments.tsv", "harmonise_log.tsv", "mvmr_input.tsv"):
        if os.path.exists(os.path.join(PROC, f)):
            z.write(os.path.join(PROC, f), f"data/proc/{f}")
# the token must never be packaged
tp = os.path.join(ROOT, ".opengwas_token")
tok = open(tp).read().strip() if os.path.exists(tp) else None
with zipfile.ZipFile(zp) as z:
    names = z.namelist()
    assert not any("opengwas_token" in n for n in names)
    if tok:
        assert not any(tok.encode() in z.read(n) for n in names), "token found in package"
if tok:
    for f in os.listdir(SUB):
        assert tok.encode() not in open(os.path.join(SUB, f), "rb").read(), f"token found in {f}"
for f in sorted(os.listdir(SUB)):
    print(f"{f:52s} {os.path.getsize(os.path.join(SUB, f)) / 1e6:6.2f} MB")
print("code zip:", os.path.relpath(zp, ROOT), "| entries:", len(names), "| token check passed")
