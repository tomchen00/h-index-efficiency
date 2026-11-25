# ======================================================================
# winsorized_hill_plot.R
#
# Purpose  : Supplementary Figures S3 (single mis-recorded count) and S4 (winsorized trade-off).
#            Draws only; the simulation lives in R/winsorized_hill_sim.R, which saves its results to data/supp/winsorized_hill.rds and then
#            sources this file. Rerun this file alone to restyle the figure in seconds.
# Produces : figures/single_outlier_influence.pdf, figures/tradeoff_winsor.pdf
# Reads    : data/supp/winsorized_hill.rds
# Requires : ggplot2 >= 4.0, ragg; R/theme_manuscript.R (xdvir and TeX for typeset labels)
# Run      : from the repository root -- Rscript R/winsorized_hill_plot.R
# ======================================================================
d <- readRDS("data/supp/winsorized_hill.rds"); df <- d$df; so <- d$so
figdir <- "figures"
source("R/theme_manuscript.R")             # house style for the paper's figures; loads ggplot2
k <- ms_col
lab <- function(tex, plain) if (ms_latex) tex else plain
# The winsorized family is ordered by q, so it takes a sequential purple ramp (dark = heavier cap); the uncapped
# plug-in is red and the empirical index ink; shapes double-encode identity.
lev <- c("raw plug-in","winsor q=0.99","winsor q=0.975","winsor q=0.95","winsor q=0.9","empirical H")  # internal level; legend text set by labels=
df$estimator <- factor(df$estimator, levels=lev)
df$kind <- ifelse(grepl("winsor",as.character(df$estimator)),"winsorized",as.character(df$estimator))
pathdf <- df[order(df$estimator), ]
pathdf <- pathdf[pathdf$kind=="winsorized", ]
cols <- setNames(c(k[["red"]], "#B9A2CF", "#9269B8", "#6A2C91", "#3E1560", k[["ink"]]), lev)
shp  <- setNames(c(15, 17, 17, 17, 17, 19), lev)
lab_tr <- if (ms_latex) c("uncapped plug-in estimator", "winsorized, $q=0.99$", "winsorized, $q=0.975$", "winsorized, $q=0.95$",
                          "winsorized, $q=0.90$", "empirical $H_n$") else
  c("uncapped plug-in estimator", "winsorized, q = 0.99", "winsorized, q = 0.975", "winsorized, q = 0.95", "winsorized, q = 0.90",
    "empirical H_n")
p <- ggplot(df, aes(contam_rmse, clean_rmse)) +
  geom_path(data=pathdf, aes(group=1), colour=k[["rule"]], linewidth=0.5, linetype="dashed") +
  geom_point(aes(colour=estimator, shape=estimator), size=2.2) +
  scale_colour_manual(values=cols, breaks=lev, labels=lab_tr, drop=FALSE) +
  scale_shape_manual(values=shp, breaks=lev, labels=lab_tr, drop=FALSE) +
  scale_x_continuous(breaks = seq(0, 50, 10)) + scale_y_continuous(breaks = seq(1.5, 3, 0.25)) +
  expand_limits(x = c(0, 50), y = c(1.5, 3)) +
  labs(x="root mean squared error, contaminated sample", y="root mean squared error, clean sample") +
  theme_manuscript() + legend_inside(0.98, 0.98)
ms_save(file.path(figdir,"tradeoff_winsor.pdf"), p, width=6.4, height=4.0)
cat("\nWrote figure:", file.path(figdir,"tradeoff_winsor.pdf"),"\n")
so$estimator <- factor(so$estimator, levels=c("raw plug-in","winsor q=0.95","empirical H"))  # internal level; legend text set by labels=
lev2 <- levels(so$estimator)
cols2 <- setNames(c(k[["red"]], k[["purple"]], k[["ink"]]), lev2)
lty2  <- setNames(c("solid", "dashed", "solid"), lev2)
lab_so <- if (ms_latex) c("Pareto plug-in estimator", "winsorized plug-in estimator, $q=0.95$", "empirical index $H_n$") else
  c("Pareto plug-in estimator", "winsorized plug-in estimator, q = 0.95", "empirical index H_n")
note <- if (ms_latex) annotate("latex", x = 1, y = 1.03, label = "bound of Proposition 1 for $H_n$", colour = k[["note"]],
                               hjust = 0, vjust = 0) else
  annotate("text", x = 1, y = 1.03, label = "bound of Proposition 1 for H_n", colour = k[["note"]], hjust = 0, vjust = 0,
           family = ms_serif, size = ms_text_size)
ps <- ggplot(so, aes(M, mean_dev, colour=estimator, linetype=estimator)) +
  geom_hline(yintercept=1, linetype="dotted", colour=k[["muted"]], linewidth=0.45) +
  geom_line() + geom_point(size=1.3, show.legend=FALSE) + note +
  scale_x_log10(breaks=10^(0:9),
                labels=ms_pow10(0:9)) +
  scale_y_continuous(breaks = seq(-0.25, 1.25, 0.25), labels = ms_minus) + expand_limits(y = c(-0.25, 1.25)) +
  scale_colour_manual(values=cols2, labels=lab_so) +
  scale_linetype_manual(values=lty2, labels=lab_so) +
  labs(x=lab("replacement count $M$", "replacement count M"), y="mean change in the estimate") +
  theme_manuscript() + legend_inside(0.02, 0.70)
ms_save(file.path(figdir,"single_outlier_influence.pdf"), ps, width=6.4, height=3.8)
cat("\nWrote figure:", file.path(figdir,"single_outlier_influence.pdf"),"\n")
