"""Build the manuscript Word files for The Prostate from manuscript/manuscript.md.

- [@Key] / [@Key1; @Key2] citations are numbered in order of first appearance and set as superscript Arabic numerals
  (AMA: after periods and commas); the reference list is AMA style from the PubMed records cached by 15_references.py.
- Tables are built from results/tables/Table1_sources.tsv and Table2_replicated.tsv as separate Word files.
- Times New Roman 12 pt, double spacing, continuous line numbers, page numbers (journal limit: 25 double-spaced pages
  for title page, abstract, body, references and figure legends).
Outputs
  manuscript/manuscript.docx                 main document (title page ... figure legends, Supporting Information list)
  manuscript/Table 1.docx, Table 2.docx      one table per file (Word table tool)
  manuscript/manuscript_with_figures.docx    reading copy with tables and the four figures at the end
Run: python3 -I code/16_build_docx.py
"""
import sys, os, re; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
import importlib.util
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_LINE_SPACING, WD_BREAK, WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

MS = os.path.join(ROOT, "manuscript")
spec = importlib.util.spec_from_file_location("refs", os.path.join(ROOT, "code", "15_references.py"))
refs = importlib.util.module_from_spec(spec); spec.loader.exec_module(refs)
REC = refs.load()

text = open(os.path.join(MS, "manuscript.md"), encoding="utf-8").read()

# ---------------------------------------------------------------- citations: numbered superscripts
order = []
for grp in re.findall(r"\[(@[^\]]+)\]", text):
    for k in re.findall(r"@([A-Za-z0-9]+)", grp):
        if k not in refs.KEYS:
            raise SystemExit(f"unknown citation key {k}")
        if k not in order:
            order.append(k)


def compress(nums):
    nums = sorted(set(nums)); out = []; i = 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(f"{nums[i]}-{nums[j]}" if j - i >= 2 else ",".join(str(n) for n in nums[i:j + 1]))
        i = j + 1
    return ",".join(out)


def cite(m):
    nums = compress([order.index(k) + 1 for k in re.findall(r"@([A-Za-z0-9]+)", m.group(1))])
    return (m.group(2) or "") + "⟦" + nums + "⟧"  # AMA: superscript follows a period or comma


text = re.sub(r"\s*\[(@[^\]]+)\]([.,])?", cite, text)
reflist = [f"{i + 1}. {refs.ama(REC[refs.KEYS[k]]).replace('--', '—')}" for i, k in enumerate(order)]

# ---------------------------------------------------------------- tables
t1 = pd.read_csv(os.path.join(TAB, "Table1_sources.tsv"), sep="\t").astype(str)
t1 = t1.apply(lambda c: c.str.replace(r"(\d)-(\d)", r"\1–\2", regex=True).str.replace("et al.", "et al", regex=False)
              .str.replace("ng/ml", "ng/mL", regex=False))
t2 = pd.read_csv(os.path.join(TAB, "Table2_replicated.tsv"), sep="\t", dtype=str)
t2 = t2[["Exposure (unit)", "Original grade", "Reanalyzed OR (95% CI)", "PSA ratio per unit (95% CI)",
         "Detection fraction (95% interval)", "Corrected OR (95% CI)", "P (corrected)", "Regrade"]]
t2 = t2.rename(columns={"Detection fraction (95% interval)": "Detection fraction (95% MC interval)"})
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
t2["P (corrected)"] = t2["P (corrected)"].str.replace(
    r"^(\d\.\d)e(-?)0*(\d+)$", lambda m: f"{m.group(1)} × 10{(m.group(2) + m.group(3)).translate(SUP)}", regex=True)
for c in t2.columns:
    t2[c] = t2[c].str.replace(r"(\d)-(\d)", r"\1–\2", regex=True).str.replace(" to -", " to −").str.replace("(-", "(−").str.replace(r"^-(?=\d)", "−", regex=True)
TABLES = [
    ("Table 1", "Data sources", t1, False,
     "BMI, body mass index; GWAS, genome-wide association study; MVP, Million Veteran Program; NHANES, National Health "
     "and Nutrition Examination Survey; OR, odds ratio; PSA, prostate-specific antigen; SD, standard deviation."),
    ("Table 2", "Detection decomposition for the 16 replicated claims (θ from cis-KLK instruments; outcome, Wang et al "
     "European meta-analysis)", t2, True,
     "ORs and PSA ratios are per unit of the exposure scale shown in parentheses. PSA ratio, ratio of geometric mean PSA "
     "per unit of exposure (exp b_PSA; BMI and smoking adjustment of the PSA GWAS undone). Detection fraction, expected "
     "detection-driven log OR (b_PSA × θ) divided by the reanalyzed log OR, with a 95% Monte Carlo (MC) interval. "
     "Corrected OR, exp(b_PCa − b_PSA × θ). Original grade as in Cui et al. Rows are ordered by detection fraction. "
     "ALT, alanine aminotransferase; BMI, body mass index; CCL2, C-C motif chemokine ligand 2; CI, confidence interval; "
     "GWAS, genome-wide association study; IGF-I, insulin-like growth factor I; LDL, low-density lipoprotein; LTL, "
     "leukocyte telomere length; OR, odds ratio; PSA, prostate-specific antigen; SD, standard deviation; θ, log OR of "
     "diagnosis per unit ln PSA."),
]


# ---------------------------------------------------------------- document helpers
def base(doc, numbered=True):
    st = doc.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(12)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    pf = st.paragraph_format; pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE; pf.space_after = Pt(0)
    for lvl, size in ((1, 13), (2, 12)):
        h = doc.styles[f"Heading {lvl}"]; h.font.name = "Times New Roman"; h.font.size = Pt(size); h.font.bold = True
        h.font.italic = lvl == 2; h.font.color.rgb = RGBColor(0, 0, 0)
        rf = h.element.rPr.find(qn("w:rFonts"))
        for a in ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme"):
            if rf is not None and rf.get(qn(a)) is not None:
                del rf.attrib[qn(a)]
        h.paragraph_format.space_before = Pt(12); h.paragraph_format.space_after = Pt(0)
        h.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = sec.top_margin = sec.bottom_margin = Cm(2.5)
    if numbered:
        line_numbers(sec)
    page_numbers(sec)


def line_numbers(sec):
    ln = OxmlElement("w:lnNumType"); ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous")
    ln.set(qn("w:distance"), "283")
    sec._sectPr.append(ln)


def page_numbers(sec):
    p = sec.footer.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run()
    for tag, txt in (("begin", None), (None, "PAGE"), ("end", None)):
        if tag:
            e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), tag)
        else:
            e = OxmlElement("w:instrText"); e.set(qn("xml:space"), "preserve"); e.text = txt
        r._r.append(e)


def add_runs(par, s, size=None):
    for part in re.split(r"(\*\*[^*]+\*\*|⟦[^⟧]+⟧)", s):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            run = par.add_run(part[2:-2]); run.bold = True
        elif part.startswith("⟦"):
            run = par.add_run(part[1:-1]); run.font.superscript = True
        else:
            run = par.add_run(part)
        if size:
            run.font.size = Pt(size)


def no_numbers(p):
    pPr = p._p.get_or_add_pPr(); pPr.append(OxmlElement("w:suppressLineNumbers"))


def set_landscape(sec):
    sec.orientation = WD_ORIENT.LANDSCAPE; sec.page_width, sec.page_height = Cm(29.7), Cm(21.0)


def table(doc, label, title, df, note):
    p = doc.add_paragraph(); add_runs(p, f"**{label}.** {title}", 11)
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE; p.paragraph_format.space_after = Pt(4)
    tb = doc.add_table(rows=len(df) + 1, cols=len(df.columns)); tb.style = "Table Grid"
    for j, c in enumerate(df.columns):
        cell = tb.cell(0, j); cell.text = ""; r = cell.paragraphs[0].add_run(str(c)); r.bold = True; r.font.size = Pt(8)
    for i, row in enumerate(df.itertuples(index=False), start=1):
        for j, v in enumerate(row):
            cell = tb.cell(i, j); cell.text = ""
            r = cell.paragraphs[0].add_run("" if pd.isna(v) or v == "nan" else str(v)); r.font.size = Pt(8)
    for row in tb.rows:
        for cell in row.cells:
            for par in cell.paragraphs:
                par.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE; no_numbers(par)
    p = doc.add_paragraph(); add_runs(p, note, 9); p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
    p.paragraph_format.space_before = Pt(4); p.paragraph_format.space_after = Pt(12)


def body(doc):
    for ln in text.split("\n"):
        if ln.startswith("# "):
            p = doc.add_paragraph(); r = p.add_run(ln[2:]); r.bold = True; r.font.size = Pt(14)
            p.paragraph_format.space_after = Pt(12)
        elif ln.startswith("## "):
            h = ln[3:]
            if h in ("Abstract", "Introduction", "References", "Figure legends"):
                doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)
            doc.add_heading(h, level=1)
        elif ln.startswith("### "):
            doc.add_heading(ln[4:], level=2)
        elif ln.strip() == "[REFERENCES]":
            for r in reflist:
                p = doc.add_paragraph(); p.add_run(r)
        elif ln.startswith("- "):
            p = doc.add_paragraph(style="List Bullet"); add_runs(p, ln[2:])
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
        elif ln.startswith(("b_adj =", "b_total =")):
            p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER; r = p.add_run(ln); r.italic = True
        elif ln.strip():
            p = doc.add_paragraph(); add_runs(p, ln)


def build_main(reading_copy):
    doc = Document(); base(doc); body(doc)
    if reading_copy:
        for label, title, df, landscape, note in TABLES:
            sec = doc.add_section(WD_SECTION.NEW_PAGE)
            if landscape:
                set_landscape(sec)
            table(doc, label, title, df, note)
        sec = doc.add_section(WD_SECTION.NEW_PAGE)
        sec.orientation = WD_ORIENT.PORTRAIT; sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
        for k, f in enumerate(["Fig1_concept_flow", "Fig2_theta", "Fig3_detection_fraction", "Fig4_regrade"], start=1):
            if k > 1:
                doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)
            doc.add_picture(os.path.join(FIG, f + ".png"), width=Cm(16))
            p = doc.add_paragraph(); add_runs(p, f"**Figure {k}**")
    out = os.path.join(MS, "manuscript_with_figures.docx" if reading_copy else "manuscript.docx")
    doc.save(out)
    return out


def build_tables():
    outs = []
    for label, title, df, landscape, note in TABLES:
        doc = Document(); base(doc, numbered=False)
        if landscape:
            set_landscape(doc.sections[0])
        table(doc, label, title, df, note)
        out = os.path.join(MS, f"{label}.docx"); doc.save(out); outs.append(out)
    return outs


if __name__ == "__main__":
    for o in [build_main(False), build_main(True)] + build_tables():
        print("wrote", os.path.relpath(o, ROOT), f"({os.path.getsize(o) / 1e6:.2f} MB)")
    src = open(os.path.join(MS, "manuscript.md"), encoding="utf-8").read()
    main = src.split("## Introduction")[1].split("## Acknowledgments")[0]
    abstract = src.split("## Abstract")[1].split("## Introduction")[0]
    wc = lambda s: len(re.findall(r"[A-Za-z0-9θρμ][A-Za-z0-9'’.\-–/]*", re.sub(r"\[@[^\]]+\]|\*\*[^*]+:\*\*", "", s)))
    print("words: abstract", wc(abstract), "| main text (Introduction to Conclusions)",
          wc(re.sub(r"^#+ .*$", "", main, flags=re.M)), "| references", len(reflist))
    run = re.search(r"\*\*Running title:\*\* (.+)", src).group(1)
    print("running title", len(run), "characters")
