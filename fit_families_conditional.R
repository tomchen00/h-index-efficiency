# ======================================================================
# fit_families_conditional.R
#
# Purpose  :
#   Tail-family fits, step 2 (conditional pipeline used for Tables 3-4). Fits
#   Pareto, zeta, Lomax, Burr XII, and lognormal to each author's cited papers
#   (x >= 1) by maximum likelihood with the likelihood left-truncated at one
#   citation; reports a two-sided tie-aware Kolmogorov-Smirnov distance against
#   the fitted conditional CDF and the model-implied h-index from the
#   conditional crossing m*S(h)/S(1-) = h. The zeta tail sum uses the full
#   Euler-Maclaurin remainder. Burr fits are bounded multi-start maximizations;
#   when the free optimum lies on a flat ridge (no interior maximum) the script
#   reports it and falls back to the tail-identified fit only if that fit is
#   within two log-likelihood units. Then runs the conditional Burr bootstrap
#   (B = 4000): parametric resamples of the cited-paper counts, refit, and the
#   variance ratio against the empirical index. Prints the earlier
#   unconditional caches alongside for comparison.
#
# Produces : data/family_fits_conditional.rds, data/burr_boot_conditional.rds; console tables
# Reads    : data/authors_data.rds, data/family_fits.rds, data/burr_boot.rds
# Requires : actuar  (CRAN package)
# Run      : from the supplement root --  Rscript fit_families_conditional.R
# Seed     : set.seed(1729)
# ======================================================================

## ensure output folders exist (some zip extractors drop empty dirs)
for (d_ in c("data","figures")) if (!dir.exists(d_)) dir.create(d_, showWarnings=FALSE)

suppressMessages({library(actuar)})
set.seed(1729)
DD  <- "data"
OUT <- "data"
A <- readRDS(file.path(DD,"authors_data.rds"))
FFold <- readRDS(file.path(DD,"family_fits.rds"))
BBold <- readRDS(file.path(DD,"burr_boot.rds"))

expected <- c(Economics="Acemo", Physics="Witten", Mathematics="Tao",
              CS_ML="Hinton", Statistics="Tibshirani", Genomics="Lander")
for (fld in names(expected)) {
  stopifnot(grepl(expected[[fld]], A[[fld]]$info$name))
  cat(sprintf("%-12s %-24s OpenAlex id=%s\n", fld, A[[fld]]$info$name, A[[fld]]$info$id))
}

hindex <- function(x){ xs<-sort(x,decreasing=TRUE); w<-which(xs>=seq_along(xs)); if(length(w)) max(w) else 0 }

# ---- two-sided, tie-aware descriptive KS against a conditional fitted CDF ----
# Fcond: vectorized conditional CDF F(x | X>=1) evaluated at the unique data values.
ks2 <- function(xc, Fcond){
  u <- sort(unique(xc)); m <- length(xc)
  Femp  <- cumsum(tabulate(match(xc, u)))/m       # Femp(u_j)
  FempL <- c(0, Femp[-length(Femp)])              # Femp(u_j^-)
  Fm <- Fcond(u)
  max(pmax(abs(Femp - Fm), abs(FempL - Fm)))
}

# ---- zeta with full Euler-Maclaurin tail ----
zeta_full <- function(s, K=2e6){ k<-1:K; sum(k^(-s)) + K^(1-s)/(s-1) - 0.5*K^(-s) }
# S^-(k) = P(X>=k) with full tail: cumulative head + EM remainder beyond Kz
make_Szeta <- function(s, hmax, Kz=2e6){
  Z <- zeta_full(s)
  kk <- 1:min(hmax+2, Kz); ck <- c(0, cumsum(kk^(-s)))     # sum_{j<k} j^-s
  em_tail <- function(a) a^(1-s)/(s-1) + 0.5*a^(-s) + (s/12)*a^(-s-1)  # sum_{j>=a} j^-s (a>=1)
  function(h){                                             # h integer-ish >=1
    h <- pmax(round(h),1)
    ifelse(h <= length(kk),
           (Z - ck[h]) / Z,                                # exact head: sum_{j>=h} = Z - sum_{j<h}
           em_tail(h) / Z)
  }
}
fit_zeta <- function(xc){
  slx <- sum(log(xc)); m <- length(xc)
  nll <- function(s) s*slx + m*log(zeta_full(s))
  s <- optimize(nll, c(1.01, 8))$minimum
  list(s=s, alpha=s-1)
}

# ---- left-truncated ML fits (support includes (0,1) => truncate at 1) ----
# family spec: d(x, par), p(x, par) (lower tail). par on log scale for positivity.
# Bounded truncated-ML fit. The truncated 3-parameter Burr has a likelihood
# ridge (shape1 -> large, shape2 -> small, scale compensating) on which the
# implied tail index alpha=shape1*shape2 degenerates while the fit barely
# improves; box bounds + multi-start + an interior/degenerate comparison make
# that visible instead of silently returning a corner solution.
tfit <- function(xc, dfun, pfun, starts_log, lower_log, upper_log, m=length(xc), fast=FALSE){
  nll <- function(lp){
    p <- exp(lp)
    S1 <- suppressWarnings(1 - pfun(1, p))
    if (!is.finite(S1) || S1 <= 1e-12) return(1e10)
    ll <- suppressWarnings(sum(dfun(xc, p, log=TRUE))) - m*log(S1)
    if (!is.finite(ll)) return(1e10)
    -ll
  }
  if (fast) starts_log <- starts_log[1]
  best <- NULL
  for (st in starts_log){
    o <- tryCatch(optim(st, nll, method="L-BFGS-B", lower=lower_log, upper=upper_log,
                        control=list(maxit=if(fast) 100 else 500)),
                  error=function(e) NULL)
    if (is.null(o)) next
    if (is.null(best) || o$value < best$value) best <- o
  }
  if (is.null(best)) return(NULL)
  at_bound <- any(abs(best$par - lower_log) < 1e-6) || any(abs(best$par - upper_log) < 1e-6)
  list(par=exp(best$par), nll=best$value, conv=best$convergence, at_bound=at_bound)
}
spec <- list(
  Lomax = list(d=function(x,p,log=FALSE) dpareto(x, shape=p[1], scale=p[2], log=log),
               p=function(x,p) ppareto(x, shape=p[1], scale=p[2]),
               starts=function(xc) list(log(c(1, stats::median(xc))), log(c(0.3, 5)), log(c(2, stats::median(xc)))),
               lower=log(c(0.02, 1e-3)), upper=log(c(50, 1e5)), alpha=function(p) p[1]),
  Burr  = list(d=function(x,p,log=FALSE) dburr(x, shape1=p[1], shape2=p[2], scale=p[3], log=log),
               p=function(x,p) pburr(x, shape1=p[1], shape2=p[2], scale=p[3]),
               starts=function(xc) list(log(c(1, 1, stats::median(xc))),
                                        log(c(0.5, 2, stats::median(xc))),
                                        log(c(2, 0.5, stats::median(xc))),
                                        log(c(length(xc)/sum(log(xc)), 1, 1))),
               lower=log(c(0.02, 0.02, 1e-3)), upper=log(c(50, 50, 1e5)), alpha=function(p) p[1]*p[2]),
  lnorm = list(d=function(x,p,log=FALSE) dlnorm(x, meanlog=log(p[1]), sdlog=p[2], log=log),
               p=function(x,p) plnorm(x, meanlog=log(p[1]), sdlog=p[2]),
               starts=function(xc) list(c(mean(log(xc)), log(sd(log(xc)))), c(0, log(1.5))),
               lower=c(log(1e-3), log(0.05)), upper=c(log(1e6), log(10)), alpha=function(p) NA))
# conditional crossing: m * S(h)/S(1-) = h  (S(1-)=1-pfun(1) for continuous)
cross_cond <- function(pfun, p, m){
  S1 <- 1 - pfun(1, p)
  f <- function(h) m * (1 - pfun(h, p))/S1 - h
  tryCatch(uniroot(f, c(1, max(2,m)))$root, error=function(e) NA)
}

res <- data.frame()
fits <- list()
for (fld in names(A)){
  x <- A[[fld]]$works$cites; x <- x[!is.na(x)]; H <- hindex(x)
  xc <- as.numeric(x[x>=1]); m <- length(xc); xmax <- max(xc)

  # Pareto (unchanged estimand; two-sided KS now)
  aP <- m/sum(log(xc))
  ksP <- ks2(xc, function(u) 1 - u^(-aP))
  hP <- m^(1/(1+aP))

  # zeta (full tail)
  fz <- fit_zeta(xc); Sz <- make_Szeta(fz$s, hmax=max(xmax, 4*m))
  ksZ <- ks2(xc, function(u) 1 - Sz(u+1))
  hZ  <- tryCatch(uniroot(function(h) m*Sz(h) - h, c(1, max(2,m)))$root, error=function(e) NA)

  # truncated continuous families
  fL <- tfit(xc, spec$Lomax$d, spec$Lomax$p, spec$Lomax$starts(xc), spec$Lomax$lower, spec$Lomax$upper)
  fB <- tfit(xc, spec$Burr$d,  spec$Burr$p,  spec$Burr$starts(xc),  spec$Burr$lower,  spec$Burr$upper)
  fN <- tfit(xc, spec$lnorm$d, spec$lnorm$p, spec$lnorm$starts(xc), spec$lnorm$lower, spec$lnorm$upper)
  # ridge diagnosis: refit Burr on the identified region (shapes <= 4.5 => alpha <= 20.25);
  # if the tail-identified fit is within ~2 log-lik units of the free fit, use it downstream
  # and record the gap -- the free optimum is then a flat degenerate ridge, not evidence.
  fBt <- tfit(xc, spec$Burr$d, spec$Burr$p, spec$Burr$starts(xc),
              log(c(0.02, 0.02, 1e-3)), log(c(4.5, 4.5, 1e5)))
  dnll_tight <- fBt$nll - fB$nll
  burr_degen <- (fB$par[1]*fB$par[2] > 20) || fB$at_bound
  if (burr_degen && dnll_tight < 2) {
    cat(sprintf("%-12s Burr free fit degenerate (alpha=%.3g, dnll_tight=%.3f) -> using tail-identified fit\n",
                fld, fB$par[1]*fB$par[2], dnll_tight))
    fB <- fBt
  } else if (burr_degen) {
    cat(sprintf("%-12s Burr free fit degenerate (alpha=%.3g) and tight fit worse by %.2f nll -> KEEPING free fit, tail index not interpretable\n",
                fld, fB$par[1]*fB$par[2], dnll_tight))
  }
  kcond <- function(spc, fit) { S1 <- 1 - spc$p(1, fit$par)
    ks2(xc, function(u) (spc$p(u, fit$par) - spc$p(1, fit$par))/S1) }
  ksL <- kcond(spec$Lomax, fL); ksB <- kcond(spec$Burr, fB); ksN <- kcond(spec$lnorm, fN)
  hL <- cross_cond(spec$Lomax$p, fL$par, m)
  hB <- cross_cond(spec$Burr$p,  fB$par, m)
  hN <- cross_cond(spec$lnorm$p, fN$par, m)

  fits[[fld]] <- list(burr=fB, m=m, xc=xc, H=H)
  res <- rbind(res, data.frame(field=fld, H=H, m=m,
    a_Pareto=round(aP,3), a_Lomax=round(fL$par[1],3), a_zeta=round(fz$alpha,3),
    a_Burr=round(fB$par[1]*fB$par[2],3),
    KS_Pareto=round(ksP,3), KS_zeta=round(ksZ,3), KS_Lomax=round(ksL,3),
    KS_Burr=round(ksB,3), KS_lnorm=round(ksN,3),
    h_Pareto=round(hP), h_zeta=round(hZ), h_Lomax=round(hL), h_Burr=round(hB), h_lnorm=round(hN),
    S1_Burr=round(1-spec$Burr$p(1,fB$par),3), Burr_degen=burr_degen, dnll_tight=round(dnll_tight,2)))
}
options(width=220)
cat("\n=== Conditional fits (two-sided KS; truncated ML; full zeta tail) ===\n")
print(res, row.names=FALSE)
cat("\n=== Earlier unconditional fits (family_fits.rds), for comparison ===\n")
print(FFold[,c("field","H","KS_Pareto","KS_zeta","KS_Lomax","KS_Burr","KS_lnorm",
               "h_Pareto","h_zeta","h_Lomax","h_Burr","h_lnorm")], row.names=FALSE)
cat(sprintf("\nConditional mean KS:  Par %.3f zeta %.3f Lom %.3f Burr %.3f lnorm %.3f\n",
  mean(res$KS_Pareto), mean(res$KS_zeta), mean(res$KS_Lomax), mean(res$KS_Burr), mean(res$KS_lnorm)))
cat(sprintf("Conditional mean |h-H|: Par %.1f zeta %.1f Lom %.1f Burr %.1f lnorm %.1f\n",
  mean(abs(res$h_Pareto-res$H)), mean(abs(res$h_zeta-res$H),na.rm=TRUE),
  mean(abs(res$h_Lomax-res$H),na.rm=TRUE), mean(abs(res$h_Burr-res$H),na.rm=TRUE),
  mean(abs(res$h_lnorm-res$H),na.rm=TRUE)))
saveRDS(res, file.path(OUT,"family_fits_conditional.rds"))

# ---- Burr bootstrap, corrected (truncated fit + conditional crossing) ----
B <- 4000
out <- data.frame()
for (fld in names(A)){
  fi <- fits[[fld]]; xc <- fi$xc; m <- fi$m; H0 <- fi$H
  x <- A[[fld]]$works$cites; x <- x[!is.na(x)]; zeros <- x[x<1]
  st_log <- log(fi$burr$par)
  # bootstrap on the SAME region the full-data fit used (tail-identified if it was chosen)
  bt_upper <- if (fi$burr$par[1] <= 4.5 && fi$burr$par[2] <= 4.5) log(c(4.5, 4.5, 1e5)) else spec$Burr$upper
  bt_lower <- spec$Burr$lower
  bH<-bP<-bB<-numeric(B); fail <- 0
  for (b in 1:B){
    xb <- sample(xc,m,replace=TRUE)
    bH[b] <- hindex(c(xb,zeros))
    aPb <- m/sum(log(xb)); bP[b] <- m^(1/(1+aPb))
    fb <- tryCatch(tfit(xb, spec$Burr$d, spec$Burr$p, list(st_log), bt_lower, bt_upper, fast=TRUE),
                   error=function(e) NULL)
    if (is.null(fb) || !is.finite(fb$nll)) {                          # cold restart once (full fit)
      fb <- tryCatch(tfit(xb, spec$Burr$d, spec$Burr$p, spec$Burr$starts(xb), bt_lower, bt_upper),
                     error=function(e) NULL)
    }
    bB[b] <- if (!is.null(fb)) cross_cond(spec$Burr$p, fb$par, m) else NA
    if (is.na(bB[b])) fail <- fail + 1
  }
  vH<-var(bH); vP<-var(bP); vB<-var(bB,na.rm=TRUE)
  hB0 <- cross_cond(spec$Burr$p, fi$burr$par, m)
  out <- rbind(out, data.frame(field=fld, H=H0, h_Burr=round(hB0,1),
    bias_Burr=round(hB0-H0,1), VarH=round(vH,1), VarBurr=round(vB,2),
    ratio_Burr=round(vB/vH,3), ratio_Pareto=round(vP/vH,3),
    conv=round(mean(!is.na(bB)),3), nfail=fail))
  cat(sprintf("%-12s H=%3d h_Burr=%6.1f bias=%5.1f | ratio Burr=%.3f Pareto=%.3f (conv %.1f%%, fail %d)\n",
              fld,H0,hB0,hB0-H0,vB/vH,vP/vH,100*mean(!is.na(bB)),fail))
}
saveRDS(out, file.path(OUT,"burr_boot_conditional.rds"))
cat("\n=== Earlier unconditional Burr bootstrap (burr_boot.rds), for comparison ===\n"); print(BBold, row.names=FALSE)
cat(sprintf("\nConditional means: |bias_Burr|=%.1f ratio_Burr=%.3f ratio_Pareto=%.3f\n",
    mean(abs(out$bias_Burr)), mean(out$ratio_Burr), mean(out$ratio_Pareto)))
cat("DONE.\n")

