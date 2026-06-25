# ======================================================================
# fit_families_boot.R
#
# Purpose  :
#   Burr bootstrap, step 3. Bootstraps the Burr XII plug-in for its calibration
#   gap and variance ratio relative to the empirical h-index (B = 4000; this is
#   the slow one).
#
# Produces : data/burr_boot.rds; console output
# Reads    : data/authors_data.rds
# Requires : fitdistrplus, actuar  (CRAN packages)
# Run      : from the supplement root --  Rscript fit_families_boot.R
# Seed     : set.seed(1729) (all stochastic scripts use this seed)
# ======================================================================

## ensure output folders exist (some zip extractors drop empty dirs)
for (d_ in c("data","figures")) if (!dir.exists(d_)) dir.create(d_, showWarnings=FALSE)

# Bootstrap the Burr-plug-in h-index (flexible RV family) vs empirical H and Pareto plug-in.
# Warm-start Burr fits from the full-data MLE for speed/stability. B modest.
suppressMessages({library(fitdistrplus); library(actuar)})
set.seed(1729)
DD <- "data"
A <- readRDS(file.path(DD,"authors_data.rds"))
hindex <- function(x){ xs<-sort(x,decreasing=TRUE); w<-which(xs>=seq_along(xs)); if(length(w)) max(w) else 0 }
cross <- function(Sfun,m) tryCatch(uniroot(function(h) m*Sfun(h)-h, c(1,max(2,m)))$root, error=function(e) NA)
burrS <- function(h,p) (1+(h/p["scale"])^p["shape2"])^(-p["shape1"])

B <- 4000
out <- data.frame()
for (fld in names(A)){
  x <- A[[fld]]$works$cites; x <- x[!is.na(x)]; H0 <- hindex(x)
  xc <- as.numeric(x[x>=1]); m <- length(xc); zeros <- x[x<1]
  f0 <- try(fitdist(xc,"burr",start=list(shape1=1,shape2=1,scale=stats::median(xc))), silent=TRUE)
  if(!inherits(f0,"fitdist")){ cat(fld,"full Burr fit failed\n"); next }
  st <- as.list(f0$estimate)
  bH<-bP<-bB<-numeric(B)
  for (b in 1:B){
    xb <- sample(xc,m,replace=TRUE)
    bH[b] <- hindex(c(xb,zeros))
    aP <- m/sum(log(xb)); bP[b] <- cross(function(h) h^(-aP), m)
    fb <- try(fitdist(xb,"burr",start=st), silent=TRUE)
    bB[b] <- if(inherits(fb,"fitdist")) cross(function(h) burrS(h,fb$estimate), m) else NA
  }
  vH<-var(bH); vP<-var(bP); vB<-var(bB,na.rm=TRUE)
  hBurr0 <- cross(function(h) burrS(h,f0$estimate), m)
  out <- rbind(out, data.frame(field=fld, H=H0, h_Burr=round(hBurr0,1),
    bias_Burr=round(hBurr0-H0,1), VarH=round(vH,1), VarBurr=round(vB,2),
    ratio_Burr=round(vB/vH,3), ratio_Pareto=round(vP/vH,3),
    conv=round(mean(!is.na(bB)),2)))
  cat(sprintf("%-12s H=%3d h_Burr=%6.1f bias=%5.1f | Var ratio Burr=%.3f Pareto=%.3f (conv %.0f%%)\n",
              fld,H0,hBurr0,hBurr0-H0,vB/vH,vP/vH,100*mean(!is.na(bB))))
}
saveRDS(out, file.path(DD,"burr_boot.rds"))
cat("\nMean: |bias_Burr|=",round(mean(abs(out$bias_Burr)),1),
    " ratio_Burr=",round(mean(out$ratio_Burr),3)," ratio_Pareto=",round(mean(out$ratio_Pareto),3),"\n")
cat("DONE.\n")
