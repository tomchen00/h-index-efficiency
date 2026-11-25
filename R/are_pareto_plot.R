# ======================================================================
# are_pareto_plot.R
#
# Purpose  : Supplementary Figure S1, Var(plug-in)/Var(H_n) for Pareto citation counts.
#            Draws only; the simulation lives in R/are_pareto_fig.R, which saves its results to data/supp/are_pareto.rds and then
#            sources this file. Rerun this file alone to restyle the figure in seconds.
# Produces : figures/are_pareto.pdf
# Reads    : data/supp/are_pareto.rds
# Requires : ggplot2 >= 4.0, ragg; R/theme_manuscript.R (xdvir and TeX for typeset labels)
# Run      : from the repository root -- Rscript R/are_pareto_plot.R
# ======================================================================
d <- readRDS("data/supp/are_pareto.rds"); grid <- d$grid; ns <- d$ns; par_tab <- d$par_tab
source("R/theme_manuscript.R")             # house style for the paper's figures; loads ggplot2
k <- ms_col
lab <- function(tex, plain) if (ms_latex) tex else plain
df <- data.frame(n=ns, ratio_sim=par_tab$VarRatio)
lev <- c("simulation", "first-order delta method", "leading term of Theorem 1")
p <- ggplot(mapping = aes(colour = what, linetype = what)) +
  geom_line(data = transform(grid, what = lev[3]), aes(n, ratio_th)) +
  geom_line(data = transform(grid, what = lev[2]), aes(n, ratio_ref)) +
  geom_point(data = transform(df, what = lev[1]), aes(n, ratio_sim), size = 1.8, show.legend = c(colour = TRUE, linetype = FALSE)) +
  scale_x_log10(breaks = c(100, 300, 1000, 3000, 10000), labels = c("100", "300", "1000", "3000", "10000")) +
  scale_y_continuous(breaks = seq(0.2, 0.6, 0.1)) + expand_limits(y = c(0.2, 0.6)) +
  scale_colour_manual(values = setNames(c(k[["ink"]], k[["blue"]], k[["muted"]]), lev), breaks = lev) +
  scale_linetype_manual(values = setNames(c("blank", "solid", "dashed"), lev), breaks = lev) +
  guides(colour = guide_legend(override.aes = list(shape = c(16, NA, NA), linetype = c("blank", "solid", "dashed"))),
         linetype = "none") +
  labs(x = lab("number of papers $n$", "number of papers n"),
       y = lab("$\\mathrm{Var}(\\hat h_n)/\\mathrm{Var}(H_n)$", "Var(plug-in) / Var(H_n)")) +
  theme_manuscript() + legend_inside(0.98, 0.98)
ms_save("figures/are_pareto.pdf", p, width = 6.2, height = 3.8)
cat("\nWrote figures/are_pareto.pdf\n")
