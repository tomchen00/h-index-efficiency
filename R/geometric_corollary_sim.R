# ======================================================================
# geometric_corollary_sim.R
#
# Purpose  :
#   Light-tailed (geometric) corollary, verification. Confirms that under a
#   light tail H_n concentrates on one or two integers (lattice-degenerate), so
#   the heavy-tailed mechanism of Theorem 1 does not apply.
#
# Produces : console diagnostics (no figure)
# Reads    : nothing (self-contained simulation)
# Requires : base R only
# Run      : from the supplement root --  Rscript geometric_corollary_sim.R
# Seed     : set.seed(1729) (all stochastic scripts use this seed)
# ======================================================================

# =====================================================================
# Geometric (light-tailed) corollary: the efficiency collapse is a HEAVY-tail
# phenomenon. For X ~ Geom(q) on {0,1,2,...}, S(h)=(1-q)^h=e^{-lambda h},
# lambda=-log(1-q). Light (non-RV) tail: elasticity alpha_n=lambda*h_n=W(lambda n)->inf.
# Predict: H_n is LATTICE-DEGENERATE (concentrates on 1-2 integers, Var(H_n) bounded/->0,
# continuized SD a_n=sqrt(h_n)/(1+alpha_n)->0), so NO Gaussian limit. The Theorem-1 MECHANISM (a
# diverging Var(H_n)) does not occur; the variance RATIO may still ->0, but via sub-lattice localization
# (plug-in resolves the crossing inside one lattice cell), not heavy-tail info-discard. Both H_n and the
# plug-in localize to O(1) of h_n.
# =====================================================================
set.seed(1729)
lambertW0 <- function(z){ w<-log(z+1); for(i in 1:80){ e<-exp(w); w<-w-(w*e-z)/(e*(w+1)) }; w }
hindex <- function(x){ xs<-sort(x,decreasing=TRUE); w<-which(xs>=seq_along(xs)); if(length(w)) max(w) else 0 }

run <- function(q, ns, B=4000){
  lam <- -log(1-q)
  cat(sprintf("\nGeometric q=%.2f  (lambda=%.4f)\n", q, lam))
  cat(sprintf("%6s %7s %7s %7s %9s %9s %9s %8s %8s\n",
              "n","h_n","alpha_n","a_n","Var(H_n)","Var(hhat)","ratio","P(mode)","#supp95"))
  for (n in ns){
    hn  <- lambertW0(lam*n)/lam            # population crossing n e^{-lam h}=h
    an  <- sqrt(hn)/(1+lam*hn)             # continuized SD scale (degeneracy diagnostic)
    H <- numeric(B); Hh <- numeric(B)
    for (b in 1:B){
      x <- rgeom(n, q)
      H[b] <- hindex(x)
      C <- sum(x); lh <- log(1+n/C)        # MLE: lambda_hat = log(1+n/C)
      Hh[b] <- lambertW0(lh*n)/lh
    }
    tb <- sort(table(H), decreasing=TRUE)
    pmode <- tb[1]/B
    # smallest # of integer values covering 95% of the H_n mass (lattice concentration)
    nsupp <- which(cumsum(tb)/B >= 0.95)[1]
    cat(sprintf("%6d %7.2f %7.2f %7.3f %9.4f %9.4f %9.3f %8.2f %8d\n",
                n, hn, lam*hn, an, var(H), var(Hh), var(Hh)/var(H), pmode, nsupp))
  }
}
run(0.3, c(100,1000,10000,100000))
run(0.1, c(100,1000,10000,100000))
cat("\nContrast: for Pareto (heavy tail) Var(H_n)->inf, ratio->0 via DIVERGING-variance collapse (Theorem 1).\n")
cat("Here Var(H_n) is tight (subsequentially 1-2 point, ->0 away from half-integers), H_n on ~1-3 integers.\n")
cat("The ratio may still ->0, but via sub-lattice/model-smoothing localization, NOT the diverging-Var(H_n) mechanism.\n")
cat("DONE.\n")
