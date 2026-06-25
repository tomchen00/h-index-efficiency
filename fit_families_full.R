# ======================================================================
# fit_families_full.R
#
# Purpose  :
#   Tail-family fits, step 2. Fits Pareto, zeta, Lomax, Burr XII, and lognormal
#   to each author's citations and reports the Kolmogorov-Smirnov distance and
#   the model-implied h-index from the crossing m*S(h)=h.
#
# Produces : data/family_fits.rds; console tables
# Reads    : data/authors_data.rds
# Requires : fitdistrplus, actuar  (CRAN packages)
# Run      : from the supplement root --  Rscript fit_families_full.R
# ======================================================================

## ensure output folders exist (some zip extractors drop empty dirs)
for (d_ in c("data","figures")) if (!dir.exists(d_)) dir.create(d_, showWarnings=FALSE)

# =====================================================================
# Full family comparison for the OpenAlex application (item 11).
# Fit, per author, the regularly varying families the efficiency theorem covers:
#   - Pareto   S(x)=x^{-alpha}                 RV, rigid single-alpha power law
#   - Lomax    S(x)=(1+x/scale)^{-alpha}       RV, flexible body, tail alpha=shape
#   - zeta     P(X=k)=k^{-s}/zeta(s)           RV (discrete), tail alpha=s-1  [Corollary 1 member]
#   - Burr XII S(x)=(1+(x/scale)^c1)^{-c2}     RV, body+curvature, tail alpha=c1*c2
# plus lognormal as the NOT-regularly-varying cautionary contrast.
# Report GoF (KS, lower=better), tail index alpha, and the model h-crossing m*S(h)=h vs realized H.
# Pareto and zeta are both rigid single-alpha power laws and fit the BODY poorly; the flexible RV
# families (Lomax, Burr) fit well -- the full landscape, not Pareto alone.
# =====================================================================
suppressMessages({library(fitdistrplus); library(actuar)})
DD <- "data"
A  <- readRDS(file.path(DD, "authors_data.rds"))

hindex <- function(x){ xs<-sort(x,decreasing=TRUE); w<-which(xs>=seq_along(xs)); if(length(w)) max(w) else 0 }
cross  <- function(Sfun, m){ f <- function(h) m*Sfun(h)-h
  tryCatch(uniroot(f, c(1, max(2,m)))$root, error=function(e) NA) }

# Riemann zeta and its log-derivative via truncation + Euler-Maclaurin tail.
zeta_fun <- function(s, K=2e6){ k<-1:K; sum(k^(-s)) + K^(1-s)/(s-1) - 0.5*K^(-s) }
# fit zeta by ML: minimize  s*sum(log x) + m*log(zeta(s))  over s>1.
fit_zeta <- function(xc){
  slx <- sum(log(xc)); m <- length(xc)
  nll <- function(s) s*slx + m*log(zeta_fun(s))
  s <- optimize(nll, c(1.01, 8))$minimum
  Z <- zeta_fun(s)
  list(s=s, alpha=s-1, Z=Z,
       S=function(h){ k<-ceiling(h); sapply(k, function(kk){
         (sum((kk:(kk+2e6))^(-s)) ) / Z }) })
}
# discrete survival S_zeta(h)=P(X>=h) computed once on a grid (cheap, h up to data max).
Szeta_grid <- function(s, Z, hmax){
  kk <- 1:(hmax+1); tailsum <- rev(cumsum(rev(kk^(-s)))) / Z   # tailsum[h]=sum_{k>=h} k^-s / Z
  function(h){ h<-pmin(pmax(round(h),1), hmax+1); tailsum[h] }
}

res <- data.frame()
for (fld in names(A)){
  x <- A[[fld]]$works$cites; x <- x[!is.na(x)]; H <- hindex(x)
  xc <- as.numeric(x[x>=1]); m <- length(xc); xmax <- max(xc)

  # --- Pareto baseline (exact, x>=1) ---
  aP <- m/sum(log(xc)); SP <- function(h) h^(-aP)
  ksP <- { z<-sort(xc); max(abs(seq_along(z)/m - (1 - z^(-aP)))) }
  hP <- cross(SP, m)

  # --- Lomax (actuar pareto: shape=alpha, scale) ---
  fL <- try(fitdist(xc,"pareto",start=list(shape=1,scale=stats::median(xc))), silent=TRUE)

  # --- zeta (discrete power law) ---
  fz <- fit_zeta(xc); Sz <- Szeta_grid(fz$s, fz$Z, xmax)
  ksZ <- { z<-sort(xc); Femp<-seq_along(z)/m; Ffit<-1 - Sz(z+1); max(abs(Femp - Ffit)) }
  hZ  <- cross(function(h) Sz(h), m)

  # --- Burr XII (shape1, shape2, scale); alpha=shape1*shape2 ---
  fB <- try(fitdist(xc,"burr",start=list(shape1=1,shape2=1,scale=stats::median(xc))), silent=TRUE)
  # --- lognormal (NOT RV) ---
  fN <- try(fitdist(xc,"lnorm"), silent=TRUE)

  ks  <- function(o) if(inherits(o,"fitdist")) tryCatch(gofstat(o)$ks, error=function(e) NA) else NA
  aL  <- if(inherits(fL,"fitdist")) unname(fL$estimate["shape"]) else NA
  aB  <- if(inherits(fB,"fitdist")) unname(fB$estimate["shape1"]*fB$estimate["shape2"]) else NA
  hL  <- if(inherits(fL,"fitdist")) cross(function(h)(1+h/fL$estimate["scale"])^(-fL$estimate["shape"]), m) else NA
  hB  <- if(inherits(fB,"fitdist")) cross(function(h)(1+(h/fB$estimate["scale"])^fB$estimate["shape2"])^(-fB$estimate["shape1"]), m) else NA
  hN  <- if(inherits(fN,"fitdist")) cross(function(h) 1-plnorm(h,fN$estimate["meanlog"],fN$estimate["sdlog"]), m) else NA

  res <- rbind(res, data.frame(field=fld, H=H, m=m,
     a_Pareto=round(aP,3), a_Lomax=round(aL,3), a_zeta=round(fz$alpha,3), a_Burr=round(aB,3),
     KS_Pareto=round(ksP,3), KS_Lomax=round(ks(fL),3), KS_zeta=round(ksZ,3),
     KS_Burr=round(ks(fB),3), KS_lnorm=round(ks(fN),3),
     h_Pareto=round(hP), h_Lomax=round(hL), h_zeta=round(hZ), h_Burr=round(hB), h_lnorm=round(hN)))
}
options(width=220)
cat("\n=== GoF (KS; lower=better) and tail index alpha ===\n")
print(res[,c("field","H","m","a_Pareto","a_Lomax","a_zeta","a_Burr",
             "KS_Pareto","KS_Lomax","KS_zeta","KS_Burr","KS_lnorm")], row.names=FALSE)
cat("\n=== Model h-crossing m*S(h)=h vs realized H ===\n")
print(res[,c("field","H","h_Pareto","h_Lomax","h_zeta","h_Burr","h_lnorm")], row.names=FALSE)
saveRDS(res, file.path(DD,"family_fits.rds"))
cat(sprintf("\nMean KS:  Pareto %.3f  Lomax %.3f  zeta %.3f  Burr %.3f  lnorm %.3f\n",
  mean(res$KS_Pareto), mean(res$KS_Lomax,na.rm=TRUE), mean(res$KS_zeta),
  mean(res$KS_Burr,na.rm=TRUE), mean(res$KS_lnorm,na.rm=TRUE)))
cat(sprintf("Mean |h_family - H|:  Pareto %.1f  Lomax %.1f  zeta %.1f  Burr %.1f  lnorm %.1f\n",
  mean(abs(res$h_Pareto-res$H)), mean(abs(res$h_Lomax-res$H),na.rm=TRUE),
  mean(abs(res$h_zeta-res$H)), mean(abs(res$h_Burr-res$H),na.rm=TRUE),
  mean(abs(res$h_lnorm-res$H),na.rm=TRUE)))
cat("\nDONE.\n")
