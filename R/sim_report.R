# ======================================================================
# sim_report.R
#
# Purpose  : Figures and tables of the simulation study and the ASN application for the manuscript.
# Reads    : data/sim/results.rds, data/sim/decisions_short_n.csv, data/sim/calibration.csv,
#            data/asn/derived/asn_summary.rds
# Produces : figures/sim_mse_ratio.png (short setting, tau = 0.3),
#            figures/sim_mse_ratio_supp.png (other settings), data/sim/report_numbers.txt
# Run      : from the repository root -- Rscript R/sim_report.R
# ======================================================================
suppressPackageStartupMessages(library(ggplot2))
r <- readRDS("data/sim/results.rds")
M <- r$mse[!r$mse$unbalanced, ]
lab_est <- c(Pareto = "Pareto", LN = "lognormal", NB2 = "NB2", Burr = "Burr XII", LNpool = "pooled lognormal",
             Burrpool = "pooled Burr XII", BurrPP = "partially pooled Burr XII", HW4 = "Hill–Weissman")
lab_fam <- c(C1 = "Burr XII control", C2 = "NB2 control", T1 = "Poisson–lognormal", T2 = "hit papers",
             T3 = "zero-inflated", T4 = "paper age", T5 = "slow tail", T7 = "dependent papers")
show <- c("Pareto", "LN", "NB2", "Burr", "LNpool", "Burrpool", "BurrPP", "HW4")
plot_mse <- function(D, file, w, h) {
  D <- D[D$estimator %in% show, ]
  D$Estimator <- factor(lab_est[D$estimator], levels = lab_est[show])
  D$Truth <- factor(lab_fam[D$family], levels = lab_fam)
  D$ratio_plot <- pmin(pmax(D$ratio, 0.25), 8)
  p <- ggplot(D, aes(n, ratio_plot, colour = Estimator, linetype = Estimator)) +
    geom_hline(yintercept = 1, colour = "grey40") + geom_line(linewidth = 0.5) + geom_point(size = 0.9) +
    scale_x_log10(breaks = c(12, 32, 128, 512)) +
    scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = c("0.25", "0.5", "1", "2", "4", "≥ 8")) +
    labs(x = expression("record length" ~ italic(n)),
         y = expression(MSE(italic(T)) / MSE(italic(H)[italic(n)]))) + theme_bw(base_size = 9) +
    theme(legend.position = "bottom", legend.title = element_blank(), panel.grid.minor = element_blank())
  if (length(unique(D$setting)) > 1 || length(unique(D$tau)) > 1)
    p <- p + facet_grid(paste0(setting, " setting, ω = ", tau) ~ Truth) else p <- p + facet_wrap(~Truth, nrow = 2)
  ggsave(file, p, width = w, height = h, dpi = 300)
}
plot_mse(M[M$setting == "short" & M$tau == 0.3, ], "figures/sim_mse_ratio.png", 7.2, 4.6)
plot_mse(M[!(M$setting == "short" & M$tau == 0.3), ], "figures/sim_mse_ratio_supp.png", 10, 7)

out <- file("data/sim/report_numbers.txt", "w")
say <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = out) }

# regime table: number of the 8 truths with ratio < 1, by estimator x n (short, tau = 0.3) and median ratio
S <- M[M$setting == "short" & M$tau == 0.3 & M$estimator %in% c(show, "Burr"), ]
for (e in c(show, "Burr")) {
  x <- S[S$estimator == e, ]
  say(sprintf("REGIME %-9s ", e), paste(sapply(c(16, 32, 64, 256, 1024), function(nn) {
    y <- x[x$n == nn, ]; sprintf("n=%d: %d/8 below, median %.2f [%.2f, %.2f]", nn, sum(y$ratio < 1), median(y$ratio), min(y$ratio), max(y$ratio)) }), collapse = " | "))
}
for (st in c("short", "long")) for (ta in c(0.1, 0.3)) for (e in c("LNpool", "Burrpool", "LN", "BurrPP", "NB2")) {
  cr <- r$crossing[r$crossing$setting == st & r$crossing$tau == ta & r$crossing$estimator == e, ]
  say(sprintf("CROSS %s tau=%.1f %-9s ", st, ta, e), paste(sprintf("%s:%s", cr$family, ifelse(is.na(cr$n_cross_up), cr$pattern, cr$n_cross_up)), collapse = "  "))
}
mx <- M[M$estimator != "H", ]
say(sprintf("MCSE max over cells (in-band ratios 0.8-1.25): %.4f", max(mx$mcse[mx$ratio > 0.8 & mx$ratio < 1.25], na.rm = TRUE)))
say(sprintf("fields per cell: min %d max %d", min(r$stopping$F), max(r$stopping$F)))
U <- r$mse[r$mse$unbalanced & r$mse$estimator %in% c("LNpool", "Burrpool", "BurrPP", "LN"), ]
for (e in unique(U$estimator)) say(sprintf("UNBAL %-9s ", e), paste(sprintf("%s/%s/n=%d:%.2f", U$family[U$estimator == e], U$setting[U$estimator == e], U$n[U$estimator == e], U$ratio[U$estimator == e])[order(U$n[U$estimator == e])], collapse = " "))
D <- r$mse[!r$mse$unbalanced & r$mse$n <= 32, ]
for (e in c(show, "Burr")) { y <- D[D$estimator == e, ]; say(sprintf("DEF %-9s n<=32 defined: median %.3f min %.3f", e, median(y$def_rate), min(y$def_rate))) }

# decision table at the observed margin (short setting, n <= 32, pooled over truths and tau)
Q <- read.csv("data/sim/decisions_short_n.csv")
Q <- Q[Q$setting == "short" & Q$margin %in% -1:1, ]
for (e in unique(Q$estimator)) for (m in -1:1) {
  y <- Q[Q$estimator == e & Q$margin == m, ]; w <- y$rows
  say(sprintf("DEC %-9s margin %+d: p_pass %.3f errH %.3f errT %.3f disagree %.3f T_right %.3f | gain range over truths x n: %.3f to %.3f",
              e, m, sum(y$p_pass * w) / sum(w), sum(y$err_H * w) / sum(w), sum(y$err_T * w) / sum(w),
              sum(y$disagree * w) / sum(w), sum(y$T_right_when_disagree * y$disagree * w, na.rm = TRUE) / sum(y$disagree * w),
              min(y$err_H - y$err_T), max(y$err_H - y$err_T)))
}
cal <- read.csv("data/sim/calibration.csv"); say("CALIB columns: ", paste(names(cal), collapse = ","))
close(out)
