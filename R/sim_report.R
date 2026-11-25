# ======================================================================
# sim_report.R
#
# Purpose  : Figure 2, Table 3, Supplementary Figure S5 and the numbers quoted in Section 5.
# Reads    : data/sim/results.rds, data/sim/decisions_short_n.csv, data/sim/calibration.csv
# Produces : figures/sim_mse_ratio.pdf (n = 16 calibration, omega = 0.3; column tau),
#            figures/sim_mse_ratio_supp.pdf (other calibrations), data/sim/report_numbers.txt
# Run      : from the repository root -- Rscript R/sim_report.R
# ======================================================================
source("R/theme_manuscript.R")             # house style for the paper's figures; loads ggplot2
r <- readRDS("data/sim/results.rds")
M <- r$mse[!r$mse$unbalanced, ]
lab_est <- c(Pareto = "Pareto", LN = "lognormal", NB2 = "NB2", Burr = "Burr XII", LNpool = "pooled lognormal",
             Burrpool = "pooled Burr XII", BurrPP = "partially pooled Burr XII", HW4 = "Hill--Weissman")
lab_fam <- c(C1 = "Burr XII benchmark", C2 = "NB2 benchmark", T1 = "Poisson--lognormal", T2 = "highly cited",
             T3 = "zero-inflated", T4 = "paper age", T5 = "slow convergence", T7 = "clustered papers")
if (!ms_latex) { lab_est <- gsub("--", "\u2013", lab_est); lab_fam <- gsub("--", "\u2013", lab_fam) }
show <- c("Pareto", "LN", "NB2", "Burr", "LNpool", "Burrpool", "BurrPP", "HW4")      # report order
show_plot <- c("Burr", "Burrpool", "BurrPP", "LN", "LNpool", "NB2", "HW4", "Pareto")  # legend order: by model family
# color = model family, line type = how the shape is obtained (solid: fitted to the record; dashed: pooled;
# dotted: partially pooled), so the eight estimators stay separable in print and in grayscale
k <- ms_col
est_col <- c(Burr = k[["purple"]], Burrpool = k[["purple"]], BurrPP = k[["purple"]], LN = k[["red"]],
             LNpool = k[["red"]], NB2 = k[["blue"]], HW4 = k[["green"]], Pareto = k[["muted"]])
est_lty <- c(Burr = "solid", Burrpool = "dashed", BurrPP = "dotted", LN = "solid", LNpool = "dashed",
             NB2 = "solid", HW4 = "solid", Pareto = "solid")
plot_mse <- function(D, file, w, h, main = TRUE) {
  D <- D[D$estimator %in% show_plot, ]
  D$Estimator <- factor(lab_est[D$estimator], levels = lab_est[show_plot])
  D$Truth <- factor(lab_fam[D$family], levels = lab_fam)
  D$ratio_plot <- pmin(pmax(D$ratio, 0.25), 8)
  p <- ggplot(D, aes(n, ratio_plot, colour = Estimator, linetype = Estimator)) +
    geom_hline(yintercept = 1, colour = k[["muted"]], linewidth = 0.4, linetype = "dashed") +
    geom_line(linewidth = 0.55) + geom_point(size = 0.7, show.legend = FALSE) +
    scale_colour_manual(values = setNames(est_col[show_plot], lab_est[show_plot])) +
    scale_linetype_manual(values = setNames(est_lty[show_plot], lab_est[show_plot])) +
    scale_x_log10(breaks = c(12, 32, 128, 1024)) +
    scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = c("0.25", "0.5", "1", "2", "4", "8")) +    # values above 8 drawn at 8 (caption)
    labs(x = if (ms_latex) "record length $n$" else "record length n",
         y = if (ms_latex) "$\\mathrm{MSE}(T)/\\mathrm{MSE}(H_n)$" else "MSE(T) / MSE(H_n)") +
    theme_manuscript()
  if (main) {
    p <- p + facet_wrap(~Truth, nrow = 3) + legend_inside(0.985, 0.04) +
      guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
  } else {
    cal_n <- ifelse(D$setting == "short", 16, 64)
    D$row <- if (ms_latex) paste0("$n=", cal_n, "$ calibration, $\\omega=", D$tau, "$") else paste0("n = ", cal_n, " calibration, \u03c9 = ", D$tau)
    p <- p %+% D + facet_grid(row ~ Truth) +
      theme(legend.position = "bottom", panel.spacing.x = unit(16, "pt"),
            strip.text.x = if (ms_latex) xdvir::element_latex(size = rel(0.9)) else element_text(size = rel(0.9)),
            strip.text.y = if (ms_latex) xdvir::element_latex(angle = -90, size = rel(0.9))
                                                       else element_text(angle = -90, size = rel(0.9))) +
      guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
  }
  ms_save(file, p, w, h)
}
plot_mse(M[M$setting == "short" & M$tau == 0.3, ], "figures/sim_mse_ratio.pdf", 7.2, 5.4)
plot_mse(M[!(M$setting == "short" & M$tau == 0.3), ], "figures/sim_mse_ratio_supp.pdf", 10, 6.2, main = FALSE)

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
lo <- mx$ratio < 0.8; hi <- mx$ratio > 1.25
say(sprintf("MCSE max over cells (ratios below 0.8): %.4f", max(mx$mcse[lo], na.rm = TRUE)))
say(sprintf("MCSE / ratio max over cells (ratios above 1.25): %.4f", max(mx$mcse[hi] / mx$ratio[hi], na.rm = TRUE)))
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
