# ======================================================================
# openalex_panel_figs.R
#
# Purpose  :
#   OpenAlex figures, step 4. Builds the two real-data figures -- the log-log
#   survival curves with fitted Pareto/Burr, and the calibration-versus-variance
#   scatter -- from the cached data and summaries (the Burr curves on the
#   log-log panel are refit here, deterministically).
#
# Produces : figures/openalex_loglog.png, figures/openalex_biasvar.png; console note
# Reads    : data/authors_data.rds, data/panel_summary.rds, data/family_fits.rds, data/burr_boot.rds
# Requires : ggplot2, tidyr, fitdistrplus, actuar  (CRAN packages)
# Run      : from the supplement root --  Rscript openalex_panel_figs.R
# ======================================================================

## ensure output folders exist (some zip extractors drop empty dirs)
for (d_ in c("data","figures")) if (!dir.exists(d_)) dir.create(d_, showWarnings=FALSE)

# Figures from the cached OpenAlex panel + family fits (no re-pull, no bootstrap).
suppressMessages({library(ggplot2); library(tidyr); library(fitdistrplus); library(actuar)})
DD <- "data"
FIG <- "figures"
A  <- readRDS(file.path(DD,"authors_data.rds"))
S  <- readRDS(file.path(DD,"panel_summary.rds"))   # Pareto plug-in bootstrap ratios
FF <- readRDS(file.path(DD,"family_fits.rds"))     # KS + h_model per family
BB <- readRDS(file.path(DD,"burr_boot.rds"))       # Burr bootstrap calibration gap + var ratio

# ---- (1) log-log survival: empirical + fitted Pareto (poor) + fitted Burr (good) ----
sd <- do.call(rbind, lapply(names(A), function(fld){
  x <- A[[fld]]$works$cites; xt <- sort(x[x>=1], decreasing=TRUE); n <- length(xt)
  xc <- as.numeric(xt); aP <- n/sum(log(xc))
  fb <- try(fitdist(xc,"burr",start=list(shape1=1,shape2=1,scale=stats::median(xc))), silent=TRUE)
  Sb <- if(inherits(fb,"fitdist")) (1+(xt/fb$estimate["scale"])^fb$estimate["shape2"])^(-fb$estimate["shape1"]) else NA
  data.frame(field=sprintf("%s (%s)", fld, sub(".* ","",A[[fld]]$info$name)),
             cites=xt, emp=(seq_len(n))/n, par=pmin(1,xt^(-aP)), bur=Sb)
}))
sl <- pivot_longer(sd, c(emp,par,bur), names_to="curve", values_to="S")
sl$curve <- factor(sl$curve, levels=c("emp","par","bur"),
                   labels=c("empirical","fitted Pareto","fitted Burr"))
p1 <- ggplot(sl, aes(cites, S, colour=curve, linetype=curve)) + geom_line(linewidth=0.7) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values=c("empirical"="#333333","fitted Pareto"="#C2453E","fitted Burr"="#2D5A3D"),name=NULL)+
  scale_linetype_manual(values=c("empirical"=1,"fitted Pareto"=2,"fitted Burr"=3),name=NULL)+
  facet_wrap(~field, scales="free") +
  labs(x="citations (log)", y="survival P(X >= x) (log)",
       title="OpenAlex citation survival with Pareto and Burr fits",
       subtitle="Empirical survival, fitted Pareto, and fitted Burr on log-log axes") +
  theme_minimal(base_size=10) + theme(legend.position="top",
                                      plot.margin=margin(6, 8, 8, 8))
ggsave(file.path(FIG,"openalex_loglog.png"), p1, width=8.4, height=5.6, dpi=150)

# ---- (2) calibration-variance: absolute relative deviation |h_model-H|/H vs bootstrap variance ratio ----
bv <- rbind(
  data.frame(field=S$field, est="empirical H", relbias=0, vratio=1),
  data.frame(field=FF$field, est="Pareto plug-in",
             relbias=abs(FF$h_Pareto-FF$H)/FF$H, vratio=S$ratio_hh[match(FF$field,S$field)]),
  data.frame(field=BB$field, est="Burr plug-in",
             relbias=abs(BB$bias_Burr)/BB$H, vratio=BB$ratio_Burr))
mn <- aggregate(cbind(relbias,vratio)~est, bv, mean)
p2 <- ggplot(bv, aes(relbias, vratio, colour=est)) +
  geom_point(alpha=0.45, size=2) +
  geom_point(data=mn, size=5, shape=18) +
  geom_hline(yintercept=1, linetype=3, colour="grey60") +
  scale_colour_manual(values=c("empirical H"="#4477AA","Pareto plug-in"="#C2453E","Burr plug-in"="#2D5A3D"),name=NULL)+
  labs(x="absolute relative deviation from realized h-index  |h_model - H| / H",
       y="bootstrap variance ratio vs empirical H",
       title="Calibration and bootstrap variance ratio on OpenAlex records",
       subtitle="Author-level points and field means for empirical, Pareto, and Burr indices") +
  theme_minimal(base_size=11) + theme(legend.position="top",
                                      plot.margin=margin(6, 8, 8, 8))
ggsave(file.path(FIG,"openalex_biasvar.png"), p2, width=7.4, height=4.8, dpi=150)

cat("Wrote openalex_loglog.png (with Burr) and openalex_biasvar.png\n")
