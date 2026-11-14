"""Claim list = the 55 significant MR associations with prostate cancer graded in Cui 2024 PLoS Med (S7 Table).

Each claim is mapped to an exposure GWAS in OpenGWAS. Rule: use the exposure GWAS of the original paper when it is in
OpenGWAS; otherwise the closest available European GWAS of the same trait (reason in `map_note`).
Instrument modes:
  tophits     genome-wide clumped hits (r2<0.001, 10 Mb) at `inst_p`; if fewer than 3 SNPs, relaxed to 5e-6 then 1e-5
  region      drug-target style: SNPs within +-100 kb of the listed genes (GRCh37, Ensembl), clumped at r2<0.1,
              analysed with an LD-correlation-aware IVW
  exclude     not re-analysed (reason in `map_note`)
Output: data/proc/claims_map.tsv
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # code/ only
from common import *
import csv, re

S7 = os.path.join(ROOT, "lit", "cui2024", "cui2024_S7_MR_associations.tsv")

# claim number (order of significant rows in S7) -> (exposure id, instrument mode, inst_p, region genes or "", note)
MAP = {
    1: ("ebi-a-GCST90029029", "tophits", 5e-8, "", "Original used Jones 2019 incl. 23andMe (n=697,828); UKB-only 'morning person' (Loh 2018) used."),
    2: ("ebi-a-GCST90000047", "tophits", 5e-8, "", "Same GWAS as original (Mills 2021)."),
    3: ("ieu-b-4877", "tophits", 5e-8, "", "GSCAN excl. 23andMe (Liu 2019); original also used UKB-based smoking GWAS."),
    4: ("ieu-a-1239", "tophits", 5e-8, "", "Lee 2018 excl. 23andMe (n=766,345); original n=1,131,881 incl. 23andMe."),
    5: ("ebi-a-GCST006099", "tophits", 5e-8, "", "Same GWAS as original (Klimentidis 2018, average acceleration)."),
    6: ("ieu-b-5106", "tophits", 5e-8, "", "Age at menarche (Day 2017, excl. 23andMe) as the puberty-timing instrument, as in the original."),
    7: ("ieu-a-89", "tophits", 5e-8, "", "Same GWAS as original (Wood 2014). Original outcome: aggressive PCa; overall PCa re-analysed."),
    8: ("ieu-b-40", "tophits", 5e-8, "", "Yengo 2018 GIANT+UKB, as in the original."),
    9: ("", "exclude", None, "", "Unfavourable-adiposity composite instrument (Martin 2021) not available as summary data."),
    10: ("met-d-Albumin", "tophits", 5e-8, "", "Same sample size as original (UKB NMR, n=115,060). Original outcome: FinnGen."),
    11: ("ieu-a-1079", "tophits", 5e-8, "", "Same GWAS as original (Evans 2013, n=2,603); instruments necessarily relaxed."),
    12: ("ebi-a-GCST90016909", "tophits", 1e-5, "", "MiBioGen, p<1e-5 as in the original."),
    13: ("ebi-a-GCST90017106", "tophits", 1e-5, "", "MiBioGen, p<1e-5 as in the original."),
    14: ("ebi-a-GCST90016960", "tophits", 1e-5, "", "MiBioGen, p<1e-5 as in the original."),
    15: ("ebi-a-GCST90016982", "tophits", 1e-5, "", "MiBioGen, p<1e-5 as in the original."),
    16: ("ebi-a-GCST90012004", "tophits", 5e-8, "", "Same GWAS as original (Folkersen 2020 SCALLOP)."),
    17: ("ebi-a-GCST90012005", "tophits", 5e-8, "", "Same GWAS as original (Folkersen 2020 SCALLOP)."),
    18: ("", "exclude", None, "", "Original outcome East Asian (Biobank Japan); outside the European-ancestry re-analysis."),
    19: ("met-d-Ala", "tophits", 5e-8, "", "Original used EPIC-Norfolk+INTERVAL (n=86,507, not in OpenGWAS); UKB NMR alanine (n=115,074) used."),
    20: ("ukb-d-30620_irnt", "tophits", 5e-8, "", "Original used Sinnott-Armstrong UKB ALT (not in OpenGWAS); Neale lab UKB ALT (same participants) used."),
    21: ("ukb-d-30770_irnt", "tophits", 5e-8, "", "Same GWAS as original (Neale lab UKB IGF-1). Original outcome: aggressive PCa."),
    22: ("prot-a-1449", "tophits", 5e-8, "", "Original used 9,331 INTERVAL participants (not in OpenGWAS); Sun 2018 (n=3,301) used."),
    23: ("prot-a-3076", "tophits", 5e-8, "", "Same GWAS as original (Sun 2018). Original outcome: UKB (Neale lab)."),
    24: ("ieu-a-1049", "tophits", 5e-8, "", "Same consortium as original (Benyamin 2014)."),
    25: ("met-d-MUFA", "tophits", 5e-8, "", "Original used Kettunen 2016 (n=13,535, not in OpenGWAS under this trait); UKB NMR MUFA used."),
    26: ("", "exclude", None, "", "Serum vitamin B12 GWAS (Grarup 2013) not in OpenGWAS and no full summary statistics found."),
    27: ("ieu-a-1051", "tophits", 5e-8, "", "Benyamin 2014 (same consortium as original GISC)."),
    28: ("ukb-d-30810_irnt", "tophits", 5e-8, "", "Original used Sinnott-Armstrong UKB phosphate (not in OpenGWAS); Neale lab UKB phosphate used."),
    29: ("met-a-388", "tophits", 5e-8, "", "Original used Long 2017 (n=1,960, not in OpenGWAS); Shin 2014 aspartate (n=7,721) used."),
    30: ("ukb-d-30780_irnt", "tophits", 5e-8, "", "Same source as original (Neale lab UKB LDL). Original outcome: early-onset PCa."),
    31: ("ukb-d-30870_irnt", "tophits", 5e-8, "", "Same source as original (Neale lab UKB TG). Original outcome: early-onset PCa."),
    32: ("ebi-a-GCST90092926", "tophits", 5e-8, "", "Original used Kettunen 2016 M-VLDL-TG; UKB NMR TG in medium VLDL (Richardson 2022) used."),
    33: ("met-c-894", "tophits", 5e-8, "", "Same GWAS as original (Kettunen 2016)."),
    34: ("met-c-863", "tophits", 5e-8, "", "Same GWAS as original (Kettunen 2016, glycoprotein acetyls)."),
    35: ("met-c-850", "tophits", 5e-8, "", "Same GWAS as original (Kettunen 2016)."),
    36: ("ebi-a-GCST004432", "tophits", 5e-8, "", "Same GWAS as original (Ahola-Olli 2017)."),
    37: ("ebi-a-GCST004421", "tophits", 5e-8, "", "Same GWAS as original (Ahola-Olli 2017)."),
    38: ("ebi-a-GCST004428", "tophits", 5e-8, "", "Same GWAS as original (Ahola-Olli 2017)."),
    39: ("ebi-a-GCST004449", "tophits", 5e-8, "", "Same GWAS as original (Ahola-Olli 2017)."),
    40: ("ebi-a-GCST004435", "tophits", 5e-8, "", "Same GWAS as original (Ahola-Olli 2017, MIG/CXCL9)."),
    41: ("ebi-a-GCST004433", "tophits", 5e-8, "", "Same GWAS as original (Ahola-Olli 2017, MIP-1b/CCL4)."),
    42: ("ebi-a-GCST004438", "tophits", 5e-8, "", "Same GWAS as original (Ahola-Olli 2017, MCP-1/CCL2)."),
    43: ("prot-a-1410", "tophits", 5e-8, "", "Same GWAS as original (Sun 2018)."),
    44: ("", "exclude", None, "", "Single-SNP MSP instrument from EPIC (n=2,254) not available as summary data."),
    45: ("met-c-920", "tophits", 5e-8, "", "Same GWAS as original (Kettunen 2016)."),
    46: ("ieu-b-4879", "tophits", 5e-8, "", "Same GWAS as original (Codd 2021)."),
    47: ("ieu-b-4868", "tophits", 5e-8, "", "Same GWAS as original (male bioavailable testosterone, UKB)."),
    48: ("ebi-a-GCST90012103", "tophits", 5e-8, "", "Free testosterone GWAS not in OpenGWAS; male bioavailable testosterone (Ruth 2020) used as proxy."),
    49: ("", "exclude", None, "", "ThyroidOmics TSH (Teumer 2018) not in OpenGWAS."),
    50: ("ukb-d-30780_irnt", "region", 5e-8, "PCSK9", "Drug target: LDL-lowering variants near PCSK9 (UKB LDL; original used GLGC 2021, not in OpenGWAS)."),
    51: ("ieu-b-38", "region", 5e-8, "CACNA1C,CACNA1D,CACNA1F,CACNA1S,CACNA2D1,CACNA2D2,CACNB1,CACNB2,CACNB3,CACNB4,CACNG1", "Drug target: SBP variants in CCB target genes (ICBP+UKB; original used UKB SBP)."),
    52: ("ukb-d-30750_irnt", "region", 5e-8, "PPARG", "Drug target: HbA1c variants near PPARG (Neale lab UKB HbA1c), as in the original."),
    53: ("ukb-b-20289", "tophits", 5e-8, "", "Same GWAS as original (UKB self-reported hyperthyroidism). Original outcome: UKB."),
    54: ("ebi-a-GCST003156", "tophits", 5e-8, "", "Same GWAS as original (Bentham 2015)."),
    55: ("ebi-a-GCST006414", "tophits", 5e-8, "", "Same GWAS as original (Nielsen 2018)."),
}

rows = list(csv.reader(open(S7), delimiter="\t"))
hdr = rows[0]
out = []
author = ""
k = 0
for r in rows[1:]:
    if r[1].strip():
        author = r[1].strip()
    grade = r[-1].strip()
    if grade == "" or grade.startswith("Null"):
        continue
    k += 1
    m = re.match(r"\s*([\d.]+)\s*\(([\d.]+)\s*,\s*([\d.]+)\)\s*([*#⸸]*)", r[9])
    or_, lo, hi, flag = (float(m.group(1)), float(m.group(2)), float(m.group(3)), m.group(4)) if m else (np.nan, np.nan, np.nan, "")
    otype = "aggressive" if "*" in flag else ("early-onset" if "#" in flag else "overall")
    if "Biobank Japan" in r[7]:
        otype = "overall (East Asian)"
    pm = re.search(r"\((\d{7,8})\)", author)
    exp_id, mode, p, genes, note = MAP[k]
    out.append(dict(claim=k, category=r[0].split("(")[0].strip().title(), source=author, source_pmid=pm.group(1) if pm else "",
                    exposure=r[2].strip(), grade=grade.split(" (")[0], orig_or=or_, orig_lo=lo, orig_hi=hi,
                    orig_unit=r[8].strip(), orig_outcome=r[7].strip(), orig_outcome_type=otype, orig_n_iv=r[4].strip(),
                    exp_id=exp_id, inst_mode=mode, inst_p=p, region_genes=genes, map_note=note))
df = pd.DataFrame(out)
assert len(df) == 55, len(df)
info = load_json(os.path.join(WORK, "opengwas_gwasinfo.json"))
df["exp_trait"] = [info[i]["trait"] if i in info else "" for i in df.exp_id]
df["exp_author"] = [f'{info[i].get("author")} {info[i].get("year")}' if i in info else "" for i in df.exp_id]
df["exp_n"] = [info[i].get("sample_size") if i in info else None for i in df.exp_id]
missing = [i for i in df.exp_id if i and i not in info]
assert not missing, missing
df.to_csv(os.path.join(PROC, "claims_map.tsv"), sep="\t", index=False)
print(df.groupby("inst_mode").size())
print(df.groupby("grade").size())
print(df[["claim", "exposure", "grade", "orig_or", "orig_outcome_type", "exp_id", "exp_trait"]].to_string(index=False))
