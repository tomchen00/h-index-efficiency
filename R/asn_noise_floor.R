# asn_noise_floor.R: decision stability of the empirical h-threshold rule in the ASN application (Section 6).
# Q: are the plug-in decision changes distinguishable from what any perturbation of H_n would produce?
# (1) noise floor: within-record nonparametric bootstrap of papers, flip rate of 1{H* >= t} vs 1{H >= t};
# (2) placebo continuous estimator H + U(-1/2, 1/2): change rate by margin;
# (3) exact binomial lower-bound rule: certify h_n >= t when P{Bin(n, t/n) >= N_n(t)} <= 0.05
#     (H_n >= t <=> N_n(t) >= t; h_n >= t <=> S^-(t) >= t/n), and the mirror rule certifying h_n < t.
# Seed: set.seed(4901). Data: prorata version (the paper's main analysis).
# Run from the repository root: Rscript R/asn_noise_floor.R
suppressPackageStartupMessages(library(jsonlite))
set.seed(4901)
A <- readRDS("data/asn/derived/asn_estimates_prorata.rds")
recs <- lapply(readLines("data/asn/derived/asn_records.jsonl"), fromJSON)
ids <- vapply(recs, `[[`, "", "id")
X <- lapply(recs, function(r) floor(as.numeric(r$x_prorata)))[match(A$id, ids)]
hemp <- function(x) { s <- sort(x, decreasing = TRUE); k <- which(s >= seq_along(s)); if (length(k)) max(k) else 0L }
stopifnot(all(vapply(X, hemp, 0) == A$H))
t <- A$S3; H <- A$H; n <- A$n; marg <- H - t
out <- file.path("data", "asn", "derived", "asn_noise_floor.out")
sink(out)
cat(sprintf("records %d; margin table:\n", nrow(A))); print(table(marg))

# (1) bootstrap flip probability per record
Bb <- 500
flip <- vapply(seq_along(X), function(u) {
  x <- X[[u]]; d0 <- H[u] >= t[u]
  mean(replicate(Bb, (hemp(sample(x, replace = TRUE)) >= t[u]) != d0))
}, 0)
cat("\n(1) within-record bootstrap of H_n (B=500): mean flip probability by margin H-t\n")
print(round(tapply(flip, marg, mean), 3))
cat(sprintf("overall expected flip rate: %.3f\n", mean(flip)))

# plug-in changes by margin (pooled Burr XII, the paper's headline estimator)
chg <- A$pass_Burrpool != A$pass_H
cat("\npooled Burr XII decision change rate by margin:\n"); print(round(tapply(chg, marg, mean), 3))
cat(sprintf("overall: %.3f\n", mean(chg)))
# do plug-in changes go where the bootstrap says H's decision is fragile?
cat("\nmean bootstrap flip prob among records the pooled Burr XII changes vs does not, within margin 0:\n")
m0 <- marg == 0
print(round(c(changed = mean(flip[m0 & chg]), unchanged = mean(flip[m0 & !chg])), 3))

# (2) placebo jitter
pl <- (H + runif(length(H), -0.5, 0.5)) >= t
cat("\n(2) placebo H + U(-1/2,1/2): change rate by margin\n"); print(round(tapply(pl != (H >= t), marg, mean), 3))
cat(sprintf("overall: %.3f; pass->fail %.3f; fail->pass %.3f\n", mean(pl != (H >= t)), mean((H >= t) & !pl), mean(!(H >= t) & pl)))
cat(sprintf("pooled Burr XII: pass->fail %.3f; fail->pass %.3f\n", mean(A$pass_H & !A$pass_Burrpool), mean(!A$pass_H & A$pass_Burrpool)))

# (3) exact binomial lower-bound rules
Nt <- vapply(seq_along(X), function(u) sum(X[[u]] >= t[u]), 0)
p0 <- pmin(t / n, 1)
pv_pass <- ifelse(t > n, 1, pbinom(Nt - 1, n, p0, lower.tail = FALSE))   # P(Bin >= Nt) under p = t/n
pv_fail <- ifelse(t > n, 0, pbinom(Nt, n, p0))                           # P(Bin <= Nt) under p = t/n
cert_pass <- pv_pass <= 0.05; cert_fail <- pv_fail <= 0.05
cat("\n(3) exact binomial rule at level 0.05: share certified pass / certified fail / undecided, by margin\n")
print(round(rbind(pass = tapply(cert_pass, marg, mean), fail = tapply(cert_fail, marg, mean),
                  undecided = tapply(!cert_pass & !cert_fail, marg, mean)), 3))
cat(sprintf("overall certified pass %.3f, certified fail %.3f, undecided %.3f\n", mean(cert_pass), mean(cert_fail), mean(!cert_pass & !cert_fail)))
cat(sprintf("records with n <= 16: undecided %.3f\n", mean((!cert_pass & !cert_fail)[n <= 16])))
sink()
cat(readLines(out), sep = "\n")

# (4) at H_n = t: upper-tail features of the applications the pooled Burr XII plug-in moves to fail vs keeps (Section 6.2).
t <- A$S3; m0 <- A$H == t
mv <- m0 & A$pass_H & !A$pass_Burrpool; st <- m0 & A$pass_Burrpool
f <- function(sel, fun) median(vapply(which(sel), function(u) fun(X[[u]], t[u], A$n[u]), 0), na.rm = TRUE)
feat <- list(
  "n (papers)" = function(x, t, n) n,
  "N(t) - t (spare papers at >= t)" = function(x, t, n) sum(x >= t) - t,
  "N(2t) (papers with >= 2t)" = function(x, t, n) sum(x >= 2 * t),
  "max count / t" = function(x, t, n) max(x) / t,
  "mean of counts >= t, / t" = function(x, t, n) mean(x[x >= t]) / t,
  "share of zero-cited papers" = function(x, t, n) mean(x == 0))
out2 <- file.path("data", "asn", "derived", "asn_tail_diag.out")
sink(out2)
cat(sprintf("margin H_n = t: %d records; pooled Burr XII moves %d to fail, keeps %d passing\n\n", sum(m0), sum(mv), sum(st)))
cat(sprintf("%-34s %10s %10s\n", "median of", "moved", "kept"))
for (nm in names(feat)) cat(sprintf("%-34s %10.2f %10.2f\n", nm, f(mv, feat[[nm]]), f(st, feat[[nm]])))
# concordance across model families at the margin
cat(sprintf("\nat H_n = t, share moved to fail: lognormal %.3f, NB2 %.3f, pooled LN %.3f, pooled Burr %.3f, PP Burr %.3f\n",
  mean(!A$pass_LN[m0]), mean(!A$pass_NB2[m0]), mean(!A$pass_LNpool[m0]), mean(!A$pass_Burrpool[m0]), mean(!A$pass_BurrPP[m0])))
both <- m0 & !A$pass_Burrpool
cat(sprintf("of those pooled Burr moves, share also moved by pooled LN %.3f, by NB2 %.3f\n", mean(!A$pass_LNpool[both]), mean(!A$pass_NB2[both])))
sink(); cat(readLines(out2), sep = "\n")
