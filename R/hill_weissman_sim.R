# ======================================================================
# hill_weissman_sim.R
#
# Purpose  :
#   Semiparametric plug-in. Compares the empirical h-index with the
#   Hill-Weissman plug-in h_k = (k X_(k+1)^a_k)^(1/(1+a_k)), where a_k is the
#   Hill estimate from the k largest counts, for k = kappa*h_n. Reports the
#   simulated variance and MSE ratios against the first-order prediction
#   (h_n/k)(1 + log^2(k/h_n)). Also reports the known-scale Pareto plug-in
#   fitted to Lomax data conditional on X >= 1, which is misspecified: its
#   pseudo-true tail index is 1/(2 log 2) rather than 1.
#
# Produces : data/hill_weissman_sim.rds; console table
# Reads    : nothing (self-contained simulation)
# Requires : base R only
# Run      : from the supplement root --  Rscript hill_weissman_sim.R
# Seed     : set.seed(1729) (all stochastic scripts use this seed)
# ======================================================================

for (d_ in c("data")) if (!dir.exists(d_)) dir.create(d_, showWarnings=FALSE)
set.seed(1729)
B <- 4000
kappas <- c(2, 4, 8, 16)

h_from_sorted <- function(xs){ w <- which(xs >= seq_along(xs)); if (length(w)) max(w) else 0 }
hw_from_sorted <- function(xs, k){                 # Hill-Weissman plug-in, xs decreasing
  X <- xs[k+1]; g <- mean(log(xs[1:k])) - log(X); a <- 1/g
  (k * X^a)^(1/(1+a))
}
pareto_cond <- function(x){                        # known-scale Pareto fitted to x >= 1 (as in Sec. 7)
  xc <- x[x >= 1]; m <- length(xc); a <- m/sum(log(xc)); m^(1/(1+a))
}

dgps <- list(
  Pareto = list(r=function(n) runif(n)^(-1),                       S=function(t) pmin(1, 1/t),       beta=Inf),
  Lomax  = list(r=function(n) runif(n)^(-1) - 1,                   S=function(t) 1/(1+t),            beta=1),
  Burr03 = list(r=function(n) (runif(n)^(-0.3) - 1)^(1/0.3),       S=function(t) (1+t^0.3)^(-10/3),  beta=0.3)
)

res <- data.frame()
for (nm in names(dgps)) for (n in c(1000, 10000)){
  d <- dgps[[nm]]
  h <- uniroot(function(t) n*d$S(t) - t, c(1e-6, n))$root
  ks <- pmin(round(kappas*h), n-1)
  H <- numeric(B); W <- matrix(NA, B, length(ks)); P <- numeric(B)
  for (b in 1:B){
    x <- d$r(n); xs <- sort(x, decreasing=TRUE)
    H[b] <- h_from_sorted(xs)
    for (j in seq_along(ks)) W[b,j] <- hw_from_sorted(xs, ks[j])
    P[b] <- pareto_cond(x)
  }
  mseH <- mean((H-h)^2)
  for (j in seq_along(ks)){
    k <- ks[j]
    res <- rbind(res, data.frame(dgp=nm, n=n, h=h, kappa=kappas[j], k=k,
      varH=var(H), mseH=mseH, varW=var(W[,j]), biasW=mean(W[,j])-h,
      var_ratio=var(W[,j])/var(H), mse_ratio=mean((W[,j]-h)^2)/mseH,
      pred=(h/k)*(1+log(k/h)^2),
      biasP=mean(P)-h, mse_ratioP=mean((P-h)^2)/mseH))
  }
  cat(sprintf("%-7s n=%6d h=%7.2f VarH=%6.2f biasH=%5.2f | Pareto(cond) bias=%7.2f MSE ratio=%8.2f\n",
              nm, n, h, var(H), mean(H)-h, mean(P)-h, mean((P-h)^2)/mseH))
}
options(width=200)
print(transform(res, h=round(h,2), varH=round(varH,2), varW=round(varW,2), biasW=round(biasW,2),
                var_ratio=round(var_ratio,3), mse_ratio=round(mse_ratio,3), pred=round(pred,3),
                biasP=round(biasP,2), mse_ratioP=round(mse_ratioP,2), mseH=round(mseH,2)), row.names=FALSE)
saveRDS(res, file.path("data", "hill_weissman_sim.rds"))
cat("DONE.\n")
