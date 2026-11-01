# ======================================================================
# asn_analysis.R
#
# Purpose  : Apply the estimators of the simulation study to the ASN 2016-18 applications that are
#            h-decisive and near the h threshold: empirical H, per-record lognormal, Burr XII and NB2
#            plug-ins, and field-pooled lognormal / Burr XII (complete pooling) and partially pooled
#            Burr XII, with donors from the same sector and level excluding the applicant's own
#            person cluster (fallback: same macro-area and level when fewer than 20 donors).
#            Records whose window contains no work are excluded. Undefined fits fall back to H.
# Reads    : data/asn/derived/asn_records.jsonl, data/asn/derived/asn_person_clusters.tsv,
#            R/sim_common.R (fitting functions)
# Produces : data/asn/derived/asn_estimates_<version>.rds  one row per application
# Run      : from the repository root -- Rscript R/asn_analysis.R [prorata|lower|upper]
# Seed     : none (deterministic)
# ======================================================================
source("R/sim_common.R")
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(TRUE)
ver <- if (length(args) >= 1) args[1] else "prorata"
min_donors <- 20L

recs <- lapply(readLines(file.path("data", "asn", "derived", "asn_records.jsonl")), fromJSON)
cl <- read.delim(file.path("data", "asn", "derived", "asn_person_clusters.tsv"), colClasses = "character")
A <- data.frame(id = vapply(recs, `[[`, "", "id"),
                sector = vapply(recs, `[[`, "", "sector"),
                level = vapply(recs, `[[`, "", "level"),
                term = vapply(recs, `[[`, "", "term"),
                I1 = vapply(recs, function(r) as.numeric(r$I1), 0), I2 = vapply(recs, function(r) as.numeric(r$I2), 0),
                I3 = vapply(recs, function(r) as.numeric(r$I3), 0), S1 = vapply(recs, function(r) as.numeric(r$S1), 0),
                S2 = vapply(recs, function(r) as.numeric(r$S2), 0), S3 = vapply(recs, function(r) as.numeric(r$S3), 0),
                stringsAsFactors = FALSE)
X <- lapply(recs, function(r) floor(as.numeric(r[[paste0("x_", ver)]])))
A$n <- vapply(X, length, 0L)
A$cluster <- cl$cluster[match(A$id, cl$id)]
A$area <- substr(A$sector, 1, 2)
keep <- A$n >= 1
if (length(args) >= 2) keep <- keep & seq_along(keep) <= as.integer(args[2])   # optional smoke-test limit
A <- A[keep, ]; X <- X[keep]
N <- nrow(A)
A$m <- vapply(X, function(x) sum(x >= 1), 0)
A$H <- vapply(X, h_emp, 0)

cat(sprintf("applications with a non-empty window: %d\n", N))
fits <- lapply(seq_len(N), function(u) {
  x <- X[[u]]
  list(LN = fit_wm(x, "LN", TRUE), Burr = fit_wm(x, "Burr", TRUE, hess = TRUE), NB2 = fit_wm(x, "NB2", FALSE))
})
defLN <- vapply(fits, function(f) f$LN$def, TRUE)
defB  <- vapply(fits, function(f) f$Burr$def, TRUE)
pdB   <- vapply(fits, function(f) isTRUE(f$Burr$pd), TRUE)
plB   <- vapply(fits, function(f) isTRUE(f$Burr$at_pl), TRUE)

na_fit <- list(th = rep(NA_real_, 3), def = FALSE, at_pl = FALSE, pd = NA)
donor_set <- function(u, ok) {
  same <- A$cluster != A$cluster[u] & ok
  d <- which(same & A$sector == A$sector[u] & A$level == A$level[u])
  src <- "sector"
  if (length(d) < min_donors) { d <- which(same & A$area == A$area[u] & A$level == A$level[u]); src <- "area" }
  list(d = d, src = src)
}
A$donor_src <- NA_character_; A$ndon_B <- NA_integer_
for (u in seq_len(N)) {
  x <- X[[u]]
  dl <- donor_set(u, defLN); db <- donor_set(u, defB); dp <- donor_set(u, defB & pdB & !plB)
  A$donor_src[u] <- db$src; A$ndon_B[u] <- length(db$d)
  fits[[u]]$LNpool <- if (length(dl$d) >= 1) {
    lsg <- median(vapply(fits[dl$d], function(f) f$LN$th[2], 0)); mu0 <- mean(log(x[x >= 1] + 0.5))
    fit_wm(x, "LN", TRUE, fixed = c(NA, lsg), starts = lapply(mu0 + c(0, -1, 1, 2), function(a) c(a, lsg)))
  } else na_fit
  fits[[u]]$Burrpool <- if (length(db$d) >= 1) {
    psi <- apply(vapply(fits[db$d], function(f) f$Burr$th[1:2], c(0, 0)), 1, median); s0 <- median(log(x[x >= 1] + 0.5))
    fit_wm(x, "Burr", TRUE, fixed = c(psi, NA), starts = lapply(s0 + c(0, -1, 1, 2), function(a) c(psi, a)))
  } else na_fit
  fits[[u]]$BurrPP <- if (length(dp$d) >= SIM$pp_min_donors) {
    pr <- pp_prior(lapply(fits[dp$d], `[[`, "Burr")); s0 <- median(log(x[x >= 1] + 0.5))
    st <- list(c(pr$mu, s0), c(pr$mu, s0 + 1), c(pr$mu, s0 - 1), c(pr$mu, s0 + 2))
    f <- fit_wm(x, "Burr", TRUE, prior = pr, starts = st, hess = TRUE); f$def <- f$def && isTRUE(f$pd); f
  } else na_fit
}

ests <- c("LN", "Burr", "NB2", "LNpool", "Burrpool", "BurrPP")
fam <- c(LN = "LN", Burr = "Burr", NB2 = "NB2", LNpool = "LN", Burrpool = "Burr", BurrPP = "Burr")
tr  <- c(LN = TRUE, Burr = TRUE, NB2 = FALSE, LNpool = TRUE, Burrpool = TRUE, BurrPP = TRUE)
for (e in ests) {
  v <- vapply(seq_len(N), function(u) {
    f <- fits[[u]][[e]]
    if (!isTRUE(f$def)) return(NA_real_)
    r <- tryCatch(plugin_wm(fam[[e]], f$th, A$n[u], A$m[u], tr[[e]]), error = function(err) NA_real_)
    if (is.finite(r)) r else NA_real_
  }, 0)
  A[[paste0("est_", e)]] <- v
  A[[paste0("def_", e)]] <- !is.na(v)
  A[[paste0("op_", e)]] <- ifelse(is.na(v), A$H, v)
  A[[paste0("pass_", e)]] <- A[[paste0("op_", e)]] >= A$S3
}
A$pass_H <- A$H >= A$S3
A$pass_official <- A$I3 >= A$S3
A$burr_tail_index <- vapply(fits, function(f) if (isTRUE(f$Burrpool$def)) exp(f$Burrpool$th[1] + f$Burrpool$th[2]) else NA_real_, 0)
A$version <- ver
out <- file.path("data", "asn", "derived", sprintf("asn_estimates_%s.rds", ver))
saveRDS(A, out)
cat(sprintf("written %s\n", out))
for (e in ests) cat(sprintf("%-9s defined %.3f  decision changes vs H: up %d  down %d\n", e, mean(A[[paste0("def_", e)]]),
                            sum(!A$pass_H & A[[paste0("pass_", e)]]), sum(A$pass_H & !A[[paste0("pass_", e)]])))
