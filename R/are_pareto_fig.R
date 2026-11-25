# ======================================================================
# are_pareto_fig.R
#
# Purpose  :
#   Supplementary Table S2 and Figure S1. Simulates the Pareto (alpha = 1)
#   variance ratio Var(plug-in)/Var(H_n) at B = 8000 replicates, computes the
#   finite-n prediction (delta-method variance over the binomial-law Var(H_n))
#   and the Theorem 1 leading rate on a dense n-grid, and plots all three.
#   The simulation block replicates the alpha = 1 block of are_simulation.R
#   draw for draw (same seed, same order), so the printed table equals
#   Supplementary Table S2 exactly.
#
# Produces : figures/are_pareto.pdf; prints Supplementary Table S2 and the finite-n check to the console
# Reads    : nothing (self-contained simulation)
# Requires : ggplot2 >= 4.0, ragg; R/theme_manuscript.R (xdvir and TeX for typeset labels)
# Run      : from the repository root -- Rscript R/are_pareto_fig.R
# Seed     : set.seed(1729)
# ======================================================================

for (d_ in c("figures")) if (!dir.exists(d_)) dir.create(d_, showWarnings=FALSE)

set.seed(1729)
hindex <- function(x){ xs<-sort(x,decreasing=TRUE); k<-which(xs>=seq_along(xs)); if(length(k)) max(k) else 0L }

## Pareto: S(x)=x^{-alpha}, x>=1.  X = U^{-1/alpha}
pareto_run <- function(n, alpha, B){
  hn <- n^(1/(1+alpha))
  H <- numeric(B); hh <- numeric(B)
  for(b in 1:B){
    X  <- runif(n)^(-1/alpha)                 # Pareto(alpha), scale 1
    H[b]  <- hindex(X)
    ahat  <- n/sum(log(X))                     # Pareto MLE = Hill estimator
    hh[b] <- n^(1/(1+ahat))
  }
  ratio_emp  <- var(hh)/var(H)
  ratio_theo <- alpha^2/(1+alpha)^2 * n^(-alpha/(1+alpha)) * log(n)^2
  data.frame(alpha=alpha, n=n, h_n=round(hn,2),
             EH=round(mean(H),2), VarH=round(var(H),3),
             Ehh=round(mean(hh),2), Varhh=round(var(hh),4),
             Bias_hh=round(mean(hh)-hn,3),
             VarRatio=signif(ratio_emp,3), VarRatio_theory=signif(ratio_theo,3),
             MSEratio=signif(mean((hh-hn)^2)/mean((H-hn)^2),3))
}

ns <- c(100,300,1000,3000,10000)
cat("=== Pareto alpha=1, B=8000 (this output is Supplementary Table S2) ===\n")
par_tab <- do.call(rbind, lapply(ns, function(n) pareto_run(n, alpha=1, B=8000)))
print(par_tab, row.names=FALSE)

## Var(H_n) from the exact binomial law P(H_n >= k) = P(Bin(n, k^-alpha) >= k).
## pbinom is evaluated only within 12 standard deviations of the mean; beyond
## that the probabilities are numerically 0/1 and far-tail evaluation is unstable.
p_ge <- function(k, n, p){
  mu <- n*p; sd <- sqrt(pmax(n*p*(1-p), 1e-12)); z <- (k-1-mu)/sd
  out <- numeric(length(k))
  out[z < -12] <- 1
  mid <- abs(z) <= 12
  out[mid] <- pbinom(k[mid]-1L, n, p[mid], lower.tail=FALSE)
  out
}
exactVarH <- function(n, alpha){
  kmax <- min(n, ceiling(8*n^(1/(1+alpha))) + 50L)
  k <- 1:kmax; pr <- p_ge(k, n, k^(-alpha))
  EH <- sum(pr); EH2 <- sum((2*k-1)*pr)
  c(EH=EH, VarH=EH2-EH^2)
}

alpha <- 1
## dense n-grid for the two deterministic curves
ngrid <- unique(round(10^seq(2, 4, by = 0.1)))
grid <- data.frame(n = ngrid)
grid$ratio_th  <- alpha^2/(1+alpha)^2 * grid$n^(-alpha/(1+alpha)) * log(grid$n)^2
grid$ratio_ref <- sapply(ngrid, function(n){
  ex <- exactVarH(n, alpha)
  dh <- n^(1/(1+alpha))*log(n)/(1+alpha)^2     # |d h_n / d alpha|, exact for Pareto
  (dh^2 * alpha^2/n) / ex["VarH"]              # first-order delta over binomial-law Var(H_n)
})
fin <- t(sapply(ns, function(n){ ex <- exactVarH(n, alpha)
  dh <- n^(1/(1+alpha))*log(n)/(1+alpha)^2
  c(n=n, exEH=unname(ex["EH"]), exVarH=unname(ex["VarH"]),
    ratio_ref=unname(dh^2*alpha^2/n/ex["VarH"])) }))
fin <- as.data.frame(fin)
cat("\n=== finite-n prediction at the table n ===\n")
print(transform(fin, exEH=round(exEH,2), exVarH=round(exVarH,3), ratio_ref=signif(ratio_ref,3)), row.names=FALSE)

## consistency: binomial-law moments track the simulated moments
stopifnot(max(abs(fin$exVarH - par_tab$VarH)/par_tab$VarH) < 0.10)
stopifnot(max(abs(fin$exEH - par_tab$EH)/par_tab$EH) < 0.02)

dir.create("data/supp", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(grid = grid, ns = ns, par_tab = par_tab), "data/supp/are_pareto.rds")   # everything the figure needs
source("R/are_pareto_plot.R")
cat("\nDONE.\n")
