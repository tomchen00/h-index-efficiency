# ======================================================================
# asn_summarize.R
#
# Purpose  : Summaries of the ASN application: for each estimator, the h-threshold decisions that
#            change relative to the empirical index (upward: H < S3 <= T; downward: H >= S3 > T), overall,
#            by level and at the margin H - S3; percentile intervals from a bootstrap that resamples
#            person clusters; agreement between the OpenAlex index and the official indicator; and
#            sensitivity to the treatment of deadline-year citations.
# Reads    : data/asn/derived/asn_estimates_{prorata,lower,upper}.rds
# Produces : data/asn/derived/asn_summary.rds and printed tables
# Run      : from the repository root -- Rscript R/asn_summarize.R
# Seed     : set.seed(20260930) for the cluster bootstrap (B = 2000)
# ======================================================================
set.seed(20260930)
B <- 2000
ests <- c("LN", "Burr", "NB2", "LNpool", "Burrpool", "BurrPP")
rd <- function(v) readRDS(file.path("data", "asn", "derived", sprintf("asn_estimates_%s.rds", v)))
A <- rd("prorata")

flips <- function(D, e) c(up = sum(!D$pass_H & D[[paste0("pass_", e)]]), down = sum(D$pass_H & !D[[paste0("pass_", e)]]))
cl <- split(seq_len(nrow(A)), A$cluster)
boot_rate <- function(D, e) {
  idx <- replicate(B, unlist(cl[sample(length(cl), replace = TRUE)], use.names = FALSE), simplify = FALSE)
  r <- vapply(idx, function(i) { f <- flips(D[i, ], e); c(f / length(i), net = (f[["up"]] - f[["down"]]) / length(i)) }, c(up = 0, down = 0, net = 0))
  apply(r, 1, quantile, c(0.025, 0.975))
}
tab <- do.call(rbind, lapply(ests, function(e) {
  f <- flips(A, e); ci <- boot_rate(A, e)
  data.frame(estimator = e, defined = mean(A[[paste0("def_", e)]]), up = f[["up"]], down = f[["down"]],
             changed_pct = 100 * sum(f) / nrow(A),
             down_pct = 100 * f[["down"]] / nrow(A), down_lo = 100 * ci[1, "down"], down_hi = 100 * ci[2, "down"],
             up_pct = 100 * f[["up"]] / nrow(A), up_lo = 100 * ci[1, "up"], up_hi = 100 * ci[2, "up"])
}))
cat("applications", nrow(A), " person clusters", length(cl), "\n")
print(tab, digits = 3, row.names = FALSE)

A$margin <- pmax(-3, pmin(3, A$H - A$S3))
bym <- do.call(rbind, lapply(c("Burrpool", "BurrPP", "LNpool"), function(e) {
  do.call(rbind, lapply(split(A, A$margin), function(D) data.frame(estimator = e, margin = D$margin[1], n = nrow(D),
    pass_H = mean(D$pass_H), up = flips(D, e)[["up"]], down = flips(D, e)[["down"]])))
}))
cat("\nflips by margin H - S3 (clipped at +/-3)\n"); print(bym, digits = 3, row.names = FALSE)

bylev <- do.call(rbind, lapply(split(A, A$level), function(D) data.frame(level = D$level[1], n = nrow(D),
  down_Burrpool = flips(D, "Burrpool")[["down"]], up_Burrpool = flips(D, "Burrpool")[["up"]],
  down_BurrPP = flips(D, "BurrPP")[["down"]], up_BurrPP = flips(D, "BurrPP")[["up"]])))
cat("\nby level (1 = full, 2 = associate)\n"); print(bylev, row.names = FALSE)

meas <- c(n = nrow(A), cor = cor(A$H, A$I3), within1 = mean(abs(A$H - A$I3) <= 1), within2 = mean(abs(A$H - A$I3) <= 2),
          med_diff = median(A$H - A$I3), agree_pass = mean(A$pass_H == A$pass_official),
          official_pass_OA_fail = sum(A$pass_official & !A$pass_H), official_fail_OA_pass = sum(!A$pass_official & A$pass_H))
cat("\nOpenAlex H vs official I3\n"); print(round(meas, 3))

sens <- do.call(rbind, lapply(c("lower", "prorata", "upper"), function(v) {
  D <- rd(v); do.call(rbind, lapply(c("Burrpool", "BurrPP"), function(e) {
    f <- flips(D, e); data.frame(version = v, estimator = e, pass_H = mean(D$pass_H), up = f[["up"]], down = f[["down"]]) }))
}))
cat("\ndeadline-year sensitivity\n"); print(sens, digits = 3, row.names = FALSE)

saveRDS(list(tab = tab, by_margin = bym, by_level = bylev, measurement = meas, sensitivity = sens),
        file.path("data", "asn", "derived", "asn_summary.rds"))
