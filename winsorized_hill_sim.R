# ======================================================================
# winsorized_hill_sim.R
#
# Purpose  :
#   Efficiency-robustness trade-off. Builds the winsorized (capped) Pareto
#   plug-in, traces its clean-vs-contaminated RMSE frontier against the
#   empirical h-index and the raw plug-in, and draws the single-outlier
#   influence comparison.
#
# Produces : figures/tradeoff_winsor.png, figures/single_outlier_influence.png; console RMSE table
# Reads    : nothing (self-contained simulation)
# Requires : ggplot2  (CRAN packages)
# Run      : from the supplement root --  Rscript winsorized_hill_sim.R
# Seed     : set.seed(1729) (all stochastic scripts use this seed)
# ======================================================================

for (d_ in c("figures")) if (!dir.exists(d_)) dir.create(d_, showWarnings=FALSE)

# Known-scale Pareto: Y=log X ~ Exp(alpha). Winsorize W=min(Y,c).
#   g(a;c)=E[min(Y,c)] = (1/a)(1-e^{-ac});  solve g(a_c;c)=mean(W) for a_c (Fisher-consistent).
#   h_c = n^{1/(1+a_c)}.  Compare to empirical H_n and raw plug-in (c=Inf).
suppressMessages(library(ggplot2))
set.seed(1729)

g_win  <- function(a, c) (1/a) * (1 - exp(-a*c))                 # E[min(Y,c)]
gp_win <- function(a, c) -(1/a^2)*(1-exp(-a*c)) + (c/a)*exp(-a*c) # d/da g
EW2    <- function(a, c) (2/a^2)*(1-exp(-a*c)) - (2*c/a)*exp(-a*c) # E[min(Y,c)^2]
V_win  <- function(a, c) EW2(a,c) - g_win(a,c)^2                  # Var(min(Y,c))
ARE_c  <- function(a, c) (a^2 * gp_win(a,c)^2) / V_win(a,c)       # ARE vs MLE (<=1)
gstar  <- function(a, c) max(g_win(a,c), c - g_win(a,c)) / abs(gp_win(a,c)) # gross-error sens.

hindex <- function(x){ xs<-sort(x,decreasing=TRUE); w<-which(xs>=seq_along(xs)); if(length(w)) max(w) else 0 }
alpha_raw <- function(x) 1/mean(log(x))                           # known-scale Pareto MLE
alpha_win <- function(x, c){                                      # solve g(a;c)=mean(min(logx,c))
  Wbar <- mean(pmin(log(x), c))
  if (Wbar <= 0 || Wbar >= c) return(alpha_raw(x))
  uniroot(function(a) g_win(a,c) - Wbar, c(1e-3, 50))$root
}
hplug <- function(a, n) n^(1/(1+a))

# ---------- clean-model efficiency table (ARE + variance ratio vs H_n) ----------
clean_table <- function(alpha, ns, qs=c(.90,.95,.975,.99), B=8000){
  cat(sprintf("\n=== Clean model, alpha=%.1f ===\n", alpha))
  cat(sprintf("ARE(c) (analytic) at log-quantile caps depends on c; reported per-n below.\n"))
  for (n in ns){
    hn <- n^(1/(1+alpha))
    H <- numeric(B); Hraw <- numeric(B); Hwin <- matrix(NA,B,length(qs))
    for (b in 1:B){
      x <- (runif(n))^(-1/alpha)                                 # Pareto(alpha), x>=1
      H[b]    <- hindex(x)
      Hraw[b] <- hplug(alpha_raw(x), n)
      cq <- quantile(log(x), qs)
      for (j in seq_along(qs)) Hwin[b,j] <- hplug(alpha_win(x, cq[j]), n)
    }
    vH <- var(H); vraw <- var(Hraw)
    vwin <- apply(Hwin,2,var)
    cat(sprintf(" n=%5d  h_n=%6.2f | Var(H)=%7.3f Var(raw)=%7.3f  ratio_raw=%.3f\n",
                n, hn, vH, vraw, vraw/vH))
    for (j in seq_along(qs))
      cat(sprintf("            q=%.3f  Var(h_c)=%7.3f  ratio_vsH=%.3f  (efficiency kept vs raw=%.2f)\n",
                  qs[j], vwin[j], vwin[j]/vH, vraw/vwin[j]))
  }
}
clean_table(1.0, c(300,1000,3000))

# ---------- trade-off: clean-RMSE vs contaminated-RMSE (alpha=1, n=1000) ----------
tradeoff <- function(alpha=1, n=1000, qs=c(.90,.95,.975,.99), eps=0.05, cfac=1e6, B=8000){
  hn <- n^(1/(1+alpha))
  est <- list(H=numeric(B), raw=numeric(B))
  win <- matrix(NA,B,length(qs)); colnames(win)<-paste0("q",qs)
  estC <- list(H=numeric(B), raw=numeric(B)); winC <- matrix(NA,B,length(qs))
  m <- ceiling(eps*n)
  for (b in 1:B){
    x <- (runif(n))^(-1/alpha)
    # clean
    est$H[b]<-hindex(x); est$raw[b]<-hplug(alpha_raw(x),n)
    cq<-quantile(log(x),qs); for(j in seq_along(qs)) win[b,j]<-hplug(alpha_win(x,cq[j]),n)
    # contaminated: inflate m random papers by factor cfac
    xc<-x; idx<-sample(n,m); xc[idx]<-xc[idx]*cfac
    estC$H[b]<-hindex(xc); estC$raw[b]<-hplug(alpha_raw(xc),n)
    cqc<-quantile(log(xc),qs); for(j in seq_along(qs)) winC[b,j]<-hplug(alpha_win(xc,cqc[j]),n)
  }
  rmse<-function(v) sqrt(mean((v-hn)^2))
  df<-data.frame(
    estimator=c("empirical H","raw plug-in",paste0("winsor q=",qs)),
    clean_rmse=c(rmse(est$H),rmse(est$raw),apply(win,2,rmse)),
    contam_rmse=c(rmse(estC$H),rmse(estC$raw),apply(winC,2,rmse)))
  df$clean_rmse_rel<-df$clean_rmse/hn; df$contam_rmse_rel<-df$contam_rmse/hn
  cat(sprintf("\n=== Trade-off (alpha=%.1f, n=%d, eps=%.2f, inflate x%g), h_n=%.1f ===\n",alpha,n,eps,cfac,hn))
  print(df[,c("estimator","clean_rmse","contam_rmse")], row.names=FALSE)
  df
}
df <- tradeoff()

figdir <- "figures"
# Winsorized family is ordered by q, so it takes a sequential blue ramp (dark = heavier cap);
# raw plug-in = vermillion; empirical H = black. Shapes double-encode identity.
lev <- c("raw plug-in","winsor q=0.99","winsor q=0.975","winsor q=0.95","winsor q=0.9","empirical H")
df$estimator <- factor(df$estimator, levels=lev)
df$kind <- ifelse(grepl("winsor",as.character(df$estimator)),"winsorized",as.character(df$estimator))
pathdf <- df[order(df$estimator), ]
pathdf <- pathdf[pathdf$kind=="winsorized", ]
cols <- c("raw plug-in"="#D55E00","winsor q=0.99"="#6BAED6","winsor q=0.975"="#4292C6",
          "winsor q=0.95"="#2171B5","winsor q=0.9"="#08306B","empirical H"="#1A1A1A")
shp  <- c("raw plug-in"=15,"winsor q=0.99"=17,"winsor q=0.975"=17,
          "winsor q=0.95"=17,"winsor q=0.9"=17,"empirical H"=19)
p <- ggplot(df, aes(contam_rmse, clean_rmse)) +
  geom_path(data=pathdf, aes(group=1), colour="grey60", linewidth=0.6, linetype=2) +
  geom_point(aes(colour=estimator, shape=estimator), size=2.4) +
  scale_colour_manual(values=cols, name=NULL, drop=FALSE) +
  scale_shape_manual(values=shp, name=NULL, drop=FALSE) +
  labs(x="RMSE under 5% inflated-outlier contamination",
       y="RMSE at clean model") +
  theme_minimal(base_size=11) +
  theme(legend.position="bottom", legend.key.height=grid::unit(1.0,"lines"),
        plot.margin=margin(6, 8, 6, 8))
ggsave(file.path(figdir,"tradeoff_winsor.png"), p, width=6.8, height=4.2, dpi=150)
cat("\nWrote figure:", file.path(figdir,"tradeoff_winsor.png"),"\n")

# ---------- single-outlier influence: one mis-recorded paper of growing magnitude ----------
# Replacing one paper's count by a contaminant of magnitude M moves the empirical H_n by at
# most one (bounded influence), drags the raw plug-in without bound (IF ~ log M), and leaves
# a winsorized plug-in bounded.
single_outlier <- function(alpha=1, n=1000, B=8000, mags=10^(0:9), qwin=0.95){
  out <- data.frame()
  for (M in mags){
    dH <- numeric(B); dRaw <- numeric(B); dWin <- numeric(B)
    for (b in 1:B){
      x   <- (runif(n))^(-1/alpha)
      H0  <- hindex(x); raw0 <- hplug(alpha_raw(x), n)
      w0  <- hplug(alpha_win(x, quantile(log(x), qwin)), n)
      xc  <- x; xc[sample(n,1)] <- M                      # one paper mis-recorded as M citations
      dH[b]   <- hindex(xc) - H0
      dRaw[b] <- hplug(alpha_raw(xc), n) - raw0
      dWin[b] <- hplug(alpha_win(xc, quantile(log(xc), qwin)), n) - w0
    }
    out <- rbind(out, data.frame(M=M,
      estimator=c("empirical H","raw plug-in","winsor q=0.95"),
      mean_dev=c(mean(dH), mean(dRaw), mean(dWin))))
  }
  out
}
so <- single_outlier()
cat(sprintf("\n=== single-outlier influence (alpha=1, n=1000): mean change in index vs contaminant M ===\n"))
print(reshape(so, idvar="M", timevar="estimator", direction="wide"), row.names=FALSE)

so$estimator <- factor(so$estimator, levels=c("raw plug-in","winsor q=0.95","empirical H"))
cols2 <- c("raw plug-in"="#D55E00","winsor q=0.95"="#0072B2","empirical H"="#1A1A1A")
lty2  <- c("raw plug-in"=1,"winsor q=0.95"=5,"empirical H"=1)
ps <- ggplot(so, aes(M, mean_dev, colour=estimator, linetype=estimator)) +
  geom_hline(yintercept=1, linetype=3, colour="grey60", linewidth=0.4) +
  geom_line(linewidth=0.7) + geom_point(size=1.6, show.legend=FALSE) +
  scale_x_log10(breaks=10^(0:9),
                labels=expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9)) +
  scale_colour_manual(values=cols2, name=NULL) +
  scale_linetype_manual(values=lty2, name=NULL) +
  annotate("text", x=10^1.2, y=1.35, label="bounded by 1", size=2.8, colour="grey40", hjust=0) +
  labs(x="citations of the single mis-recorded paper",
       y="mean change in the index") +
  theme_minimal(base_size=11) +
  theme(legend.position="top", plot.margin=margin(6, 8, 6, 8))
ggsave(file.path(figdir,"single_outlier_influence.png"), ps, width=6.6, height=4.0, dpi=150)
cat("\nWrote figure:", file.path(figdir,"single_outlier_influence.png"),"\n")

# analytic gross-error sensitivity + ARE across cap c (alpha=1)
cat("\n=== analytic ARE(c) and gross-error sensitivity (alpha=1) ===\n")
for (c in c(1,2,3,4,5,8)) cat(sprintf(" c=%2d  ARE=%.3f  gamma*=%.2f\n", c, ARE_c(1,c), gstar(1,c)))
cat("\nDONE.\n")
