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
# Produces : figures/crossing.pdf; prints the numbers to the console
# Reads    : nothing (self-contained simulation)
# Requires : ggplot2 >= 4.0, patchwork, ragg; R/theme_manuscript.R (xdvir and TeX for typeset labels)
# Run      : from the repository root -- Rscript R/crossing_fig.R
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

# ---- figure (house style: R/theme_manuscript.R) ----------------------------------------------------------------
source("R/theme_manuscript.R")
suppressPackageStartupMessages(library(patchwork))
k <- ms_col
math <- function(x, y, tex, plotmath, colour = k[["ink"]], hjust = 0, vjust = 0.5) {
  if (ms_latex) annotate("latex", x = x, y = y, label = tex, colour = colour, hjust = hjust, vjust = vjust)
  else annotate("text", x = x, y = y, label = plotmath, parse = TRUE, colour = colour, hjust = hjust, vjust = vjust,
                family = ms_serif, size = ms_text_size)
}

# (a) one record
R <- 26; xs <- sort(X, decreasing = TRUE)[1:R]
pts <- data.frame(r = 1:R, y = pmin(xs, 30))
rr <- seq(0.5, R, length.out = 400)
lev <- c("ranked citation counts", "population curve", "fitted Pareto curve", "diagonal")
cur <- rbind(data.frame(r = rr, y = (n / rr)^(1 / alpha), what = lev[2]),
             data.frame(r = rr, y = (n / rr)^(1 / ax), what = lev[3]))
diag_df <- data.frame(r = c(0, 30), y = c(0, 30), what = lev[4])
col_a <- setNames(c(k[["note"]], k[["ink"]], k[["red"]], k[["muted"]]), lev)
lty_a <- setNames(c("blank", "solid", "dashed", "dotted"), lev)
pa <- ggplot(mapping = aes(r, y, colour = what, linetype = what)) +
  geom_line(data = diag_df, linewidth = 0.45) +
  geom_line(data = cur) +
  geom_point(data = transform(pts, what = lev[1]), size = 1.1, show.legend = c(colour = TRUE, linetype = FALSE)) +
  annotate("segment", x = Hx, xend = Hx, y = 0, yend = Hx, colour = k[["note"]], linewidth = 0.35) +
  annotate("point", x = hn, y = hn, shape = 21, fill = k[["paper"]], colour = k[["ink"]], size = 2.4, stroke = 0.7) +
  annotate("point", x = hx, y = hx, shape = 24, fill = k[["paper"]], colour = k[["red"]], size = 2.1, stroke = 0.7) +
  math(Hx + 0.4, 1.4, sprintf("$H_n=%d$", Hx), sprintf("italic(H)[italic(n)]==%d", Hx), colour = k[["note"]]) +
  annotate("segment", x = hn - 0.4, y = hn - 0.6, xend = 7.4, yend = 4.4, colour = k[["muted"]], linewidth = 0.3) +
  math(7.3, 3.6, "$h_n=10$", "italic(h)[italic(n)]==10", hjust = 1) +
  annotate("segment", x = hx + 0.5, y = hx + 0.8, xend = 15.2, yend = 15.4, colour = k[["red"]], linewidth = 0.3) +
  math(15.4, 15.9, sprintf("$\\hat h_n=%.1f$", hx), sprintf("hat(italic(h))[italic(n)]==%.1f", hx),
       colour = k[["red"]]) +
  scale_colour_manual(values = col_a, breaks = lev) + scale_linetype_manual(values = lty_a, breaks = lev) +
  scale_x_continuous(breaks = seq(0, 25, 5)) + scale_y_continuous(breaks = seq(0, 30, 10)) +
  coord_cartesian(xlim = c(0, R), ylim = c(0, 30), expand = FALSE) +
  labs(title = "(a) 1 record", x = "rank $r$", y = "citations") +
  theme_manuscript() + legend_inside(0.98, 0.98) +
  guides(colour = guide_legend(override.aes = list(shape = c(16, NA, NA, NA), linetype = unname(lty_a))),
         linetype = "none")

# (b) 8000 records
hb <- hist(H, breaks = seq(min(H) - 0.5, max(H) + 0.5, by = 1), plot = FALSE)
bars <- data.frame(xmin = hb$breaks[-length(hb$breaks)] + 0.08, xmax = hb$breaks[-1] - 0.08, ymax = hb$counts / B)
d <- density(hh, adjust = 1.2); dens <- data.frame(x = d$x, y = d$y)
dens <- dens[dens$x >= 4 & dens$x <= 17, ]
ytop <- max(bars$ymax)
pb <- ggplot() +
  geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = ymax), fill = "#D9D9D9") +
  geom_line(data = dens, aes(x, y), colour = k[["red"]]) +
  annotate("segment", x = hn, xend = hn, y = 0, yend = 0.40, colour = k[["muted"]], linetype = "dotted", linewidth = 0.45) +
  math(hn + 0.25, 0.405, "$h_n=10$", "italic(h)[italic(n)]==10", vjust = 0) +
  math(4.3, 0.31, sprintf("$H_n$, variance %.2f", var(H)), sprintf("italic(H)[italic(n)]*', variance %.2f'", var(H)),
       colour = k[["note"]]) +
  math(11.9, 0.25, sprintf("$\\hat h_n$, variance %.2f", var(hh)),
       sprintf("hat(italic(h))[italic(n)]*', variance %.2f'", var(hh)), colour = k[["red"]]) +
  scale_x_continuous(breaks = seq(4, 16, 2)) +
  scale_y_continuous(breaks = seq(0, 0.4, 0.1), labels = c("0", "0.1", "0.2", "0.3", "0.4")) +
  coord_cartesian(xlim = c(4, 17), ylim = c(0, 0.44), expand = FALSE) +
  labs(title = sprintf("(b) %d records", B), x = "estimate of $h_n$", y = "density") +
  theme_manuscript()

ms_save("figures/crossing.pdf", pa + pb + plot_layout(widths = c(1, 1)), width = 7.2, height = 3.1)

cat("\nExact P(H_n >= t), Pareto alpha = 1, n = 100:\n")
for (t in 7:13) cat(sprintf("  t = %2d: %.4f\n", t, pbinom(t - 1, n, t^(-alpha), lower.tail = FALSE)))
cat(sprintf("P(7 <= H_n <= 12) = %.4f; P(H_n <= 9) = %.4f\n",
            pbinom(6, n, 7^-alpha, lower.tail = FALSE) - pbinom(12, n, 13^-alpha, lower.tail = FALSE),
            1 - pbinom(9, n, 10^-alpha, lower.tail = FALSE)))
cat("\nP{Bin(n, t/n) >= t}:\n")
for (nt in list(c(16, 7), c(100, 10)))
  cat(sprintf("  n = %3d, t = %2d: %.4f\n", nt[1], nt[2], pbinom(nt[2] - 1, nt[1], nt[2] / nt[1], lower.tail = FALSE)))

# Section 2, plug-in pass probabilities (continuous Pareto, n = 100, threshold t = 10).
# hat h_n = n^{1/(1 + hat alpha)} >= t iff hat alpha <= log n / log t - 1 iff sum log X >= n / (log n / log t - 1),
# with sum log X ~ Gamma(n, rate alpha); the empirical rule passes with P{Bin(n, t^-alpha) >= t}.
cat("\nPass probabilities at t = 10, n = 100:\n")
cat(sprintf("  h_n = 10 (alpha = 1): empirical %.4f, plug-in %.4f\n",
            pbinom(9, n, 10^-alpha, lower.tail = FALSE), pgamma(n, shape = n, rate = alpha, lower.tail = FALSE)))
for (h in c(9, 11)) {
  a <- log(n) / log(h) - 1                      # alpha giving h_n = h at n = 100
  pH <- pbinom(9, n, 10^(-a), lower.tail = FALSE)
  pP <- pgamma(n, shape = n, rate = a, lower.tail = FALSE)
  err <- if (h < 10) c(pH, pP) else c(1 - pH, 1 - pP)
  cat(sprintf("  h_n = %d (alpha = %.3f): error empirical %.3f, plug-in %.3f\n", h, a, err[1], err[2]))
}
