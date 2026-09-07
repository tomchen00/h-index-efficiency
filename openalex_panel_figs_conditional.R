# ======================================================================
# openalex_panel_figs_conditional.R
#
# Purpose  :
#   OpenAlex figures, step 3 (conditional pipeline). Builds the two real-data
#   figures -- the log-log survival curves with the fitted Pareto and Burr
#   conditional survival functions, and the calibration-versus-variance
#   scatter -- from the cached data and the conditional fit summaries (the Burr
#   curves on the log-log panel are refit here, deterministically, with the
#   same truncated likelihood as fit_families_conditional.R).
#
# Produces : figures/openalex_loglog.png, figures/openalex_biasvar.png; console note
# Reads    : data/authors_data.rds, data/panel_summary.rds,
#            data/family_fits_conditional.rds, data/burr_boot_conditional.rds
# Requires : ggplot2, tidyr, actuar  (CRAN packages)
# Run      : from the supplement root --  Rscript openalex_panel_figs_conditional.R
# Seed     : none needed (deterministic given the caches)
# ======================================================================

## ensure output folders exist (some zip extractors drop empty dirs)
for (d_ in c("data","figures")) if (!dir.exists(d_)) dir.create(d_, showWarnings=FALSE)

suppressMessages({library(ggplot2); library(tidyr); library(actuar)})
DD  <- "data"
CD  <- "data"
OUT <- "figures"
A  <- readRDS(file.path(DD,"authors_data.rds"))
S  <- readRDS(file.path(DD,"panel_summary.rds"))
FF <- readRDS(file.path(CD,"family_fits_conditional.rds"))
BB <- readRDS(file.path(CD,"burr_boot_conditional.rds"))
# --- truncated-ML Burr fit, identical to fit_families_conditional.R ---
tfit <- function(xc, dfun, pfun, starts_log, lower_log, upper_log, m=length(xc)){
  nll <- function(lp){ p <- exp(lp); S1 <- suppressWarnings(1 - pfun(1, p))
    if (!is.finite(S1) || S1 <= 1e-12) return(1e10)
    ll <- suppressWarnings(sum(dfun(xc, p, log=TRUE))) - m*log(S1); if (!is.finite(ll)) return(1e10); -ll }
  best <- NULL
  for (st in starts_log){ o <- tryCatch(optim(st, nll, method="L-BFGS-B", lower=lower_log, upper=upper_log, control=list(maxit=500)), error=function(e) NULL)
    if (is.null(o)) next; if (is.null(best) || o$value < best$value) best <- o }
  at_bound <- any(abs(best$par - lower_log) < 1e-6) || any(abs(best$par - upper_log) < 1e-6)
  list(par=exp(best$par), nll=best$value, at_bound=at_bound)
}
dB <- function(x,p,log=FALSE) dburr(x, shape1=p[1], shape2=p[2], scale=p[3], log=log)
pB <- function(x,p) pburr(x, shape1=p[1], shape2=p[2], scale=p[3])
starts <- function(xc) list(log(c(1,1,stats::median(xc))), log(c(0.5,2,stats::median(xc))), log(c(2,0.5,stats::median(xc))), log(c(length(xc)/sum(log(xc)),1,1)))
sd <- do.call(rbind, lapply(names(A), function(fld){
  x <- A[[fld]]$works$cites; x <- x[!is.na(x)]; xt <- sort(x[x>=1], decreasing=TRUE); n <- length(xt); xc <- as.numeric(xt)
  aP <- n/sum(log(xc))
  fB <- tfit(xc, dB, pB, starts(xc), log(c(0.02,0.02,1e-3)), log(c(50,50,1e5)))
  fBt <- tfit(xc, dB, pB, starts(xc), log(c(0.02,0.02,1e-3)), log(c(4.5,4.5,1e5)))
  degen <- (fB$par[1]*fB$par[2] > 20) || fB$at_bound
  if (degen && (fBt$nll - fB$nll) < 2) fB <- fBt
  cat(sprintf("%-12s alpha_Burr=%.3f (cache %.3f) degen=%s\n", fld, fB$par[1]*fB$par[2], FF$a_Burr[FF$field==fld], degen))
  S1 <- 1 - pB(1, fB$par)
  Sb <- (1 - pB(xt, fB$par))/S1            # conditional survival P(X>=x | X>=1), continuous
  data.frame(field=sprintf("%s (%s)", fld, sub(".* ","",A[[fld]]$info$name)),
             cites=xt, emp=seq_len(n)/n, par=pmin(1, xt^(-aP)), bur=pmin(1, Sb))
}))
sl <- pivot_longer(sd, c(emp,par,bur), names_to="curve", values_to="S")
sl$curve <- factor(sl$curve, levels=c("emp","par","bur"), labels=c("empirical","fitted Pareto","fitted Burr"))
p1 <- ggplot(sl, aes(cites, S, colour=curve, linetype=curve)) + geom_line(linewidth=0.6) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values=c("empirical"="#1A1A1A","fitted Pareto"="#D55E00","fitted Burr"="#0072B2"),name=NULL)+
  scale_linetype_manual(values=c("empirical"=1,"fitted Pareto"=2,"fitted Burr"=3),name=NULL)+
  facet_wrap(~field, scales="free") + labs(x="citations", y=expression(P(X>=x))) +
  theme_minimal(base_size=10) + theme(legend.position="top", plot.margin=margin(6,8,6,8))
ggsave(file.path(OUT,"openalex_loglog.png"), p1, width=8.0, height=5.2, dpi=600)
bv <- rbind(
  data.frame(field=S$field, est="empirical H", relbias=0, vratio=1),
  data.frame(field=FF$field, est="Pareto plug-in", relbias=abs(FF$h_Pareto-FF$H)/FF$H, vratio=S$ratio_hh[match(FF$field,S$field)]),
  data.frame(field=BB$field, est="Burr plug-in", relbias=abs(BB$bias_Burr)/BB$H, vratio=BB$ratio_Burr))
mn <- aggregate(cbind(relbias,vratio)~est, bv, mean); print(mn)
p2 <- ggplot(bv, aes(relbias, vratio, colour=est)) + geom_point(alpha=0.55, size=1.8) + geom_point(data=mn, size=3.6, shape=18) +
  geom_hline(yintercept=1, linetype=3, colour="grey60", linewidth=0.4) +
  scale_colour_manual(values=c("empirical H"="#1A1A1A","Pareto plug-in"="#D55E00","Burr plug-in"="#0072B2"),name=NULL)+
  labs(x="absolute relative deviation from realized h-index", y="bootstrap variance ratio vs empirical H") +
  theme_minimal(base_size=11) + theme(legend.position="top", plot.margin=margin(6,8,6,8))
ggsave(file.path(OUT,"openalex_biasvar.png"), p2, width=6.8, height=4.4, dpi=600)
cat("wrote conditional figures to", OUT, "\n")
