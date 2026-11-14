# How much of the Mendelian randomization evidence on prostate cancer reflects PSA-driven detection? A reanalysis of 50 published causal claims

Analysis code and derived tables for the manuscript submitted to *The Prostate* (2026).

Repository: https://github.com/xiaoyangguo2024/pca-mr-psa-detection

## Data sources

All inputs are publicly available summary-level or de-identified survey data. Sources are listed in `results/tables/Table1_sources.tsv` (manuscript Table 1) and include:

- **GWAS Catalog**
  - [GCST90274714](https://www.ebi.ac.uk/gwas/studies/GCST90274714) — prostate cancer, Wang et al. European meta-analysis
  - [GCST90461907](https://www.ebi.ac.uk/gwas/studies/GCST90461907) — PSA levels
- **IEU OpenGWAS** — exposure and additional outcome GWAS (identifiers in Table 1 and Table S1)
- **NHANES** — public survey files used to estimate the coefficients of ln PSA on BMI and smoking ([NHANES](https://wwwn.cdc.gov/nchs/nhanes/))

Raw GWAS files and the LD reference panel are not shipped in this repository (they are large). Download them with the scripts below.

OpenGWAS API calls require a personal token from https://api.opengwas.io. Place the token in a file named `.opengwas_token` in the project root. The token is not included here; `.gitignore` excludes it.

## Requirements

- Python 3.9+ with the packages in `requirements.txt` (`pandas`, `numpy`, `scipy`, `statsmodels`, `matplotlib`, `python-docx`, `PyMuPDF`)
- R (developed with R 4.5.3) with `data.table`, `MendelianRandomization` (0.10.0) and `jsonlite`. These are the packages used by `08_mr.R` and `12_mvmr.R`. `code/install_r_packages.R` can also install a broader set of optional MR packages into a project-local `Rlib/` directory.

Run scripts from the **project root**, for example:

```bash
python3 -I code/03_nhanes_psa_covariates.py
```

Paths in the code are relative to the project root (`code/common.py` defines `data/raw`, `data/ref`, `data/proc`, `work` and `results/`).

## Analysis steps

| Step | Script | Description |
| --- | --- | --- |
| 1 | `code/01_download_outcomes.sh`, `code/01b_download_one.sh` | Download GWAS Catalog harmonized files (MD5 checked) |
| 2 | `code/02_download_ldref.py` | Download the IEU 1000 Genomes v3 European LD reference |
| 3 | `code/03_nhanes_psa_covariates.py` | NHANES coefficients of ln PSA on BMI and smoking |
| 4 | `code/04_claims_map.py` | Map claims from Cui et al. 2024 to exposure GWAS |
| 5 | `code/05_instruments.py`, `code/06b_psa_instruments.py` | Exposure and PSA instruments (OpenGWAS clumping) |
| 6 | `code/06_lookups.py` | Look up instruments in the outcome GWAS |
| 7 | `code/07_harmonise.py` | Harmonize alleles; undo BMI/smoking adjustment of the PSA GWAS |
| 8 | `code/08_mr.R` | Univariable Mendelian randomization |
| 9 | `code/09_detection.py` | θ, detection decomposition and regrading |
| 10 | `code/10_simulation.py` | Analytic model and MR simulation |
| 11 | `code/11_mvmr_prep.py`, `code/12_mvmr.R` | Multivariable MR |
| 12 | `code/13_figures.py`, `code/14_tables.py` | Figures and tables |
| 13 | `code/15_references.py`, `code/16_build_docx.py` | PubMed reference list and manuscript Word files |
| 14 | `code/17_submission.py` | Assemble a journal submission package |

`16_build_docx.py` and `17_submission.py` generate manuscript and submission files from the analysis outputs. They do not contain the manuscript text itself. The manuscript and cover letter are not included in this repository.

Helper modules: `code/common.py` (shared paths), `code/ogwas.py` (OpenGWAS API client), `code/gwascat.py` (GWAS Catalog helpers).

## Derived tables

Shipped under `results/tables/` (analysis outputs used in the manuscript) and `data/proc/` (intermediate files needed to regenerate later steps):

**`results/tables/`**

- `Table1_sources.tsv`, `Table2_replicated.tsv` — manuscript Tables 1 and 2
- `mr_all.tsv` — univariable MR results
- `theta.tsv` — effect of PSA on diagnosis (θ)
- `detection_PSA_KLK.tsv`, `detection_PSA_gw.tsv` — detection decomposition (cis-*KLK* and genome-wide PSA instruments)
- `regrade_PSA_KLK.tsv`, `regrade_PSA_gw.tsv` — regrading of published claims
- `DF_by_estimator.tsv` — detection fractions by MR estimator
- `mvmr.tsv` — multivariable MR
- `nhanes_flow.tsv`, `nhanes_psa_covariates.tsv` — NHANES sample flow and covariate coefficients
- `sim_theta.tsv`, `sim_mr.tsv` — analytic model and simulation
- `summary_numbers.json` — key numbers used in figures and text

**`data/proc/`**

- `claims_map.tsv` — published claims mapped to exposure GWAS
- `instrument_summary.tsv`, `psa_instruments.tsv` — instrument lists and summaries
- `harmonise_log.tsv` — harmonization log
- `mvmr_input.tsv` — multivariable MR input

## License

This repository is released under the [MIT License](LICENSE). Copyright (c) 2026 Xiaoyang Guo.

## Citation

Guo X, You H, Li J, Liu Q, Sun Z, Wen B. How much of the Mendelian randomization evidence on prostate cancer reflects PSA-driven detection? A reanalysis of 50 published causal claims. Manuscript submitted to *The Prostate*. 2026. Code: https://github.com/xiaoyangguo2024/pca-mr-psa-detection
