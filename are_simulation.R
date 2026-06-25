# ======================================================================
# are_simulation.R
#
# Purpose  :
#   Theorem 1, by simulation. Compares the empirical h-index H_n with the Pareto
#   maximum-likelihood plug-in, showing the variance ratio Var(plug-in)/Var(H_n)
#   shrink toward zero as n grows; also shows the light-tailed geometric case
#   where H_n is lattice-degenerate (the boundary of the phenomenon).
#
# Produces : figures/are_pareto.png; prints the Pareto variance-ratio table to the console
# Reads    : nothing (self-contained simulation)
# Requires : base R only
# Run      : from the supplement root --  Rscript are_simulation.R
# Seed     : set.seed(1729) (all stochastic scripts use this seed)
# ======================================================================

## ensure output folders exist (some zip extractors drop empty dirs)
for (d_ in c("data","figures")) if (!dir.exists(d_)) dir.create(d_, showWarnings=FALSE)

## ARE simulation: empirical h-index H_n vs Pareto MLE plug-in h_hat = n^{1/(1+alpha_hat)}.
## Claim (Theorem 1): Var(h_hat)/Var(H_n) ~ (alpha^2/(1+alpha)^2) n^{-alpha/(1+alpha)} log^2 n -> 0.
## Also: geometric (light tail) -> H_n lattice-degenerate, Var(H_n) bounded (caveat regime).
set.seed(1729)
lambertW0 <- function(z){ w<-log(z+1); for(i in 1:80){ e<-exp(w); w<-w-(w*e-z)/(e*(w+1)) }; w }
hindex <- function(x){ xs<-sort(x,decreasing=TRUE); k<-which(xs>=seq_along(xs)); if(length(k)) max(k) else 0L }

## ---- Pareto: S(x)=x^{-alpha}, x>=1.  X = U^{-1/alpha} ----
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

cat("=== PARETO: empirical H_n vs MLE plug-in (B=8000) ===\n")
par_tab <- do.call(rbind, lapply(c(100,300,1000,3000,10000),
                                 function(n) pareto_run(n, alpha=1, B=8000)))
print(par_tab, row.names=FALSE)
cat("\n(alpha=2)\n")
par_tab2 <- do.call(rbind, lapply(c(100,300,1000,3000),
                                  function(n) pareto_run(n, alpha=2, B=8000)))
print(par_tab2, row.names=FALSE)

## ---- Geometric: S(h)=(1-q)^h -> light tail, lattice-degenerate ----
geom_run <- function(n, q, B){
  a<--log(1-q); hn<-lambertW0(a*n)/a
  H<-numeric(B); for(b in 1:B){ H[b]<-hindex(rgeom(n,q)) }
  data.frame(q=q, n=n, h_n=round(hn,2), EH=round(mean(H),2),
             VarH=round(var(H),3), modal_frac=round(max(table(H))/B,3))
}
cat("\n=== GEOMETRIC: H_n is lattice-degenerate (VarH stays bounded) ===\n")
geo_tab <- do.call(rbind, lapply(c(50,100,200,500,1000,2000,5000),
                                 function(n) geom_run(n, q=0.2, B=8000)))
print(geo_tab, row.names=FALSE)

## ---- figure: Pareto variance ratio -> 0 ----
fig <- "figures/are_pareto.png"
png(fig, width=1500, height=1100, res=200)
ns <- c(100,300,1000,3000,10000)
emp <- par_tab$VarRatio; theo <- par_tab$VarRatio_theory
plot(ns, emp, log="xy", type="b", pch=19, col="#2D5A3D", lwd=2,
     xlab="n (papers)", ylab=expression(Var(hat(h))/Var(H[n])),
     main="Pareto citations: variance ratio by sample size")
lines(ns, theo, type="b", pch=1, lty=2, col="#C2453E", lwd=2)
legend("topright", c("empirical (simulation)","theory  ~ n^{-1/2} log^2 n"),
       col=c("#2D5A3D","#C2453E"), pch=c(19,1), lty=c(1,2), lwd=2, bty="n")
dev.off()
cat("\nfigure saved:", fig, "\n")
