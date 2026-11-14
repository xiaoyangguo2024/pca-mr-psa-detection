# Two-sample MR for every exposure set x outcome in data/proc/harmonised.tsv (MendelianRandomization 0.10.0).
#   1 SNP        Wald ratio (first-order SE)
#   2-3 SNPs     IVW fixed effect
#   >=4 SNPs     IVW multiplicative random effects (residual SE floored at 1), MR-Egger, weighted median,
#                weighted mode (MBE, simple SE); Cochran's Q and the Egger intercept are reported
#   region keys  (drug-target claims 50-52 and the cis-KLK PSA instruments) IVW with the EUR LD correlation matrix (correl = TRUE) as the primary estimate
# Output: results/tables/mr_all.tsv
suppressPackageStartupMessages({library(data.table); library(MendelianRandomization); library(jsonlite)})
set.seed(20261113)
h <- fread("data/proc/harmonised.tsv")
region <- c("claim50", "claim51", "claim52", "PSA_KLK")

ld_for <- function(key, d) {
  # LD matrix from OpenGWAS ld/matrix (1000G EUR); names are rsid_A1_A2 and correlations refer to A1.
  # Returns the rows of d that are in the matrix with a matching allele, and R signed for each SNP's effect allele.
  obj <- fromJSON(sprintf("work/instruments/%s.json", key), simplifyVector = TRUE)
  if (is.null(obj$ld)) return(NULL)
  nm <- obj$ld$snplist
  M <- matrix(as.numeric(unlist(obj$ld$matrix)), nrow = length(nm), byrow = TRUE)
  sp <- strsplit(nm, "_", fixed = TRUE)
  rs <- vapply(sp, `[`, "", 1); a1 <- vapply(sp, `[`, "", 2); a2 <- vapply(sp, `[`, "", 3)
  idx <- match(d$rsid, rs)
  s <- ifelse(toupper(d$ea) == a1[idx], 1, ifelse(toupper(d$ea) == a2[idx], -1, NA))
  keep <- !is.na(idx) & !is.na(s)
  list(keep = keep, R = M[idx[keep], idx[keep], drop = FALSE] * outer(s[keep], s[keep]))
}

row <- function(key, outcome, method, n, b, se, p, extra = list()) {
  c(list(key = key, outcome = outcome, method = method, nsnp = n, b = b, se = se, p = p), extra)
}

out <- list()
for (kk in unique(h$key)) {
  for (oo in unique(h[key == kk]$outcome)) {
    d <- h[key == kk & outcome == oo]
    d <- d[is.finite(bx) & is.finite(by) & is.finite(sx) & is.finite(sy) & sy > 0]
    n <- nrow(d)
    if (n == 0) next
    if (n == 1) {
      b <- d$by / d$bx; se <- abs(d$sy / d$bx)
      out[[length(out) + 1]] <- row(kk, oo, "Wald ratio", 1, b, se, 2 * pnorm(-abs(b / se)))
      next
    }
    if (kk %in% region) {
      L <- ld_for(kk, d)
      if (!is.null(L) && sum(L$keep) >= 2) {
        dd <- d[L$keep]
        inp <- mr_input(bx = dd$bx, bxse = dd$sx, by = dd$by, byse = dd$sy, correlation = L$R, snps = dd$rsid)
        f <- tryCatch(mr_ivw(inp, correl = TRUE), error = function(e) NULL)
        if (!is.null(f)) out[[length(out) + 1]] <- row(kk, oo, "IVW (LD-correlated)", nrow(dd), f@Estimate, f@StdError,
                                                        f@Pvalue, list(Q = f@Heter.Stat[1], Q_p = f@Heter.Stat[2]))
      }
    }
    inp <- mr_input(bx = d$bx, bxse = d$sx, by = d$by, byse = d$sy, snps = d$rsid)
    f <- mr_ivw(inp, model = if (n < 4) "fixed" else "random")
    out[[length(out) + 1]] <- row(kk, oo, "IVW", n, f@Estimate, f@StdError, f@Pvalue,
                                  list(Q = f@Heter.Stat[1], Q_p = f@Heter.Stat[2]))
    if (n >= 4) {
      e <- tryCatch(mr_egger(inp), error = function(e) NULL)
      if (!is.null(e)) out[[length(out) + 1]] <- row(kk, oo, "MR-Egger", n, e@Estimate, e@StdError.Est, e@Pvalue.Est,
                                                      list(egger_int = e@Intercept, egger_int_se = e@StdError.Int,
                                                           egger_int_p = e@Pvalue.Int))
      m <- tryCatch(mr_median(inp, weighting = "weighted", iterations = 2000), error = function(e) NULL)
      if (!is.null(m)) out[[length(out) + 1]] <- row(kk, oo, "Weighted median", n, m@Estimate, m@StdError, m@Pvalue)
      w <- tryCatch(mr_mbe(inp, weighting = "weighted", stderror = "simple", iterations = 2000), error = function(e) NULL)
      if (!is.null(w)) out[[length(out) + 1]] <- row(kk, oo, "Weighted mode", n, w@Estimate, w@StdError, w@Pvalue)
    }
  }
  cat(kk, "done\n")
}
res <- rbindlist(out, fill = TRUE)
fwrite(res, "results/tables/mr_all.tsv", sep = "\t")
cat("rows", nrow(res), "\n")
