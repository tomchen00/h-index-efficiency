# ======================================================================
# crossing_fig.R
#
# Purpose  :
#   Figure 1 of the paper. Pareto citations (alpha = 1, scale 1) for a
#   record of n = 100 papers, so that h_n = 10.
#   (a) One simulated record: ranked citation counts X_(r), the population
#       curve r -> F^{-1}(1 - r/n) = (n/r)^{1/alpha}, and the fitted Pareto
#       curve (n/r)^{1/alpha_hat}. H_n, h_n and the plug-in hat h_n are the
#       points where these cross the diagonal.
#   (b) Distribution of H_n and of the plug-in over B = 8000 records.
#   The simulation uses the seed and draw order of are_pareto_fig.R, so the
#   variances printed here equal the n = 100 row of Supplementary Table S2.
#   Also prints the exact threshold probabilities quoted in Sections 1-2:
#   P(H_n >= t) = P{Bin(n, S(t)) >= t} and P{Bin(n, t/n) >= t}.
#
# Produces : figures/crossing.png; prints the numbers to the console
# Reads    : nothing (self-contained simulation)
# Requires : base R (grDevices with cairo)
# Run      : from the supplement root --  Rscript R/crossing_fig.R
# Seed     : set.seed(1729)
# ======================================================================

if (!dir.exists("figures")) dir.create("figures", showWarnings = FALSE)
set.seed(1729)
n <- 100; alpha <- 1; hn <- n^(1 / (1 + alpha)); B <- 8000
hindex <- function(x) { xs <- sort(x, decreasing = TRUE); k <- which(xs >= seq_along(xs)); if (length(k)) max(k) else 0L }

H <- numeric(B); hh <- numeric(B); first <- NULL
for (b in 1:B) {
  X <- runif(n)^(-1 / alpha)                  # Pareto(alpha), scale 1
  H[b] <- hindex(X)
  ahat <- n / sum(log(X))                     # Pareto MLE = Hill estimator
  hh[b] <- n^(1 / (1 + ahat))
  if (b == 1) first <- X
}
cat(sprintf("Var(H_n) = %.3f, Var(plug-in) = %.3f, ratio = %.3f\n", var(H), var(hh), var(hh) / var(H)))

X <- first; Hx <- hindex(X); ax <- n / sum(log(X)); hx <- n^(1 / (1 + ax))
cat(sprintf("Panel (a) record: H_n = %d, alpha_hat = %.3f, plug-in = %.2f\n", Hx, ax, hx))

col_fit <- "#B03A2E"
png("figures/crossing.png", width = 7.2, height = 3.3, units = "in", res = 300, type = "cairo")
par(mfrow = c(1, 2), mar = c(3.6, 3.8, 1.6, 0.6), mgp = c(2.3, 0.7, 0), las = 1, cex = 0.8)

R <- 26; xs <- sort(X, decreasing = TRUE)[1:R]
plot(1:R, pmin(xs, 30), pch = 16, cex = 0.6, col = "grey35", xlim = c(0, R), ylim = c(0, 30),
     xlab = expression("rank" ~ italic(r)), ylab = "citations")
title(main = "(a) one record", font.main = 1, cex.main = 1, adj = 0)
abline(0, 1, lty = 3)
rr <- seq(0.5, R, length.out = 400)
lines(rr, (n / rr)^(1 / alpha), lwd = 1.8)
lines(rr, (n / rr)^(1 / ax), lwd = 1.8, lty = 2, col = col_fit)
points(hn, hn, pch = 1, cex = 1.4, lwd = 1.6)
points(hx, hx, pch = 2, cex = 1.2, lwd = 1.6, col = col_fit)
segments(Hx, 0, Hx, Hx, col = "grey35")
text(Hx, 1.3, bquote(italic(H)[italic(n)] == .(Hx)), pos = 4, col = "grey20")
segments(hn, hn, 7, 4.2, col = "grey40"); text(5.4, 3.2, bquote(italic(h)[italic(n)] == 10), pos = 4)
segments(hx, hx, 13.5, 22, col = col_fit)
text(13.5, 22, bquote(hat(italic(h))[italic(n)] == .(sprintf("%.1f", hx))), pos = 4, col = col_fit)
legend("topright", bty = "n", cex = 0.85,
       legend = c("ranked citation counts", "population curve", "fitted Pareto curve", "diagonal"),
       lty = c(NA, 1, 2, 3), pch = c(16, NA, NA, NA), col = c("grey35", "black", col_fit, "black"),
       lwd = c(NA, 1.8, 1.8, 1))

hb <- hist(H, breaks = seq(min(H) - 0.5, max(H) + 0.5, by = 1), plot = FALSE)
plot(NA, xlim = c(4, 17), ylim = c(0, 0.42), xlab = expression("estimate of" ~ italic(h)[italic(n)]),
     ylab = "density")
title(main = sprintf("(b) %d records", B), font.main = 1, cex.main = 1, adj = 0)
rect(hb$breaks[-length(hb$breaks)] + 0.08, 0, hb$breaks[-1] - 0.08, hb$counts / B, col = "grey80", border = NA)
d <- density(hh, adjust = 1.2); lines(d$x, d$y, lwd = 1.8, col = col_fit)
abline(v = hn, lty = 3)
legend("topright", bty = "n", cex = 0.85,
       legend = c(bquote(italic(H)[italic(n)] * ", variance" ~ .(sprintf("%.2f", var(H)))),
                  bquote(hat(italic(h))[italic(n)] * ", variance" ~ .(sprintf("%.2f", var(hh)))),
                  bquote(italic(h)[italic(n)] == 10)),
       fill = c("grey80", NA, NA), border = NA, lty = c(NA, 1, 3), col = c(NA, col_fit, "black"), lwd = c(NA, 1.8, 1))
invisible(dev.off())

cat("\nExact P(H_n >= t), Pareto alpha = 1, n = 100:\n")
for (t in 7:13) cat(sprintf("  t = %2d: %.4f\n", t, pbinom(t - 1, n, t^(-alpha), lower.tail = FALSE)))
cat(sprintf("P(7 <= H_n <= 12) = %.4f; P(H_n <= 9) = %.4f\n",
            pbinom(6, n, 7^-alpha, lower.tail = FALSE) - pbinom(12, n, 13^-alpha, lower.tail = FALSE),
            1 - pbinom(9, n, 10^-alpha, lower.tail = FALSE)))
cat("\nP{Bin(n, t/n) >= t}:\n")
for (nt in list(c(16, 7), c(100, 10)))
  cat(sprintf("  n = %3d, t = %2d: %.4f\n", nt[1], nt[2], pbinom(nt[2] - 1, nt[1], nt[2] / nt[1], lower.tail = FALSE)))
