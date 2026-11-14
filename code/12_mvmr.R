# Multivariable MR: exposure X and ln PSA (BMI/smoking adjustment undone) on prostate cancer (Wang 2023).
# IVW (multiplicative random effects) and MV-Egger; conditional F-statistics are not computable without the
# exposure-exposure covariance, so instrument strength is described by the univariable F (see instrument_summary.tsv).
# Output: results/tables/mvmr.tsv (direct effect of X not through PSA, and the effect of PSA given X)
suppressPackageStartupMessages({library(data.table); library(MendelianRandomization)})
d <- fread("data/proc/mvmr_input.tsv")
out <- list()
for (k in unique(d$key)) {
  s <- d[key == k]
  if (nrow(s) < 6 || sum(s$src == "X") < 2) next
  inp <- mr_mvinput(bx = cbind(s$bX, s$bP), bxse = cbind(s$sX, s$sP), by = s$bY, byse = s$sY,
                    exposure = c("X", "lnPSA"), snps = s$rsid)
  f <- mr_mvivw(inp, model = "random")
  out[[length(out) + 1]] <- data.table(set = k, method = "MV-IVW", term = c("X", "lnPSA"), nsnp = nrow(s),
                                       n_x = sum(s$src == "X"), n_psa = sum(s$src == "PSA"),
                                       b = f@Estimate, se = f@StdError, p = f@Pvalue)
  e <- tryCatch(mr_mvegger(inp, orientate = 1), error = function(e) NULL)
  if (!is.null(e)) out[[length(out) + 1]] <- data.table(set = k, method = "MV-Egger", term = c("X", "lnPSA"),
                                                        nsnp = nrow(s), n_x = sum(s$src == "X"), n_psa = sum(s$src == "PSA"),
                                                        b = e@Estimate, se = e@StdError.Est, p = e@Pvalue.Est,
                                                        egger_int = e@Intercept, egger_int_p = e@Pvalue.Int)
}
res <- rbindlist(out, fill = TRUE)
setnames(res, "set", "key")  # data.table() reserves the argument name "key"
fwrite(res, "results/tables/mvmr.tsv", sep = "\t")
print(res[method == "MV-IVW"])
