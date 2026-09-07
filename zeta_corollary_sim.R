# ======================================================================
# zeta_corollary_sim.R
#
# Purpose  :
#   Discrete (zeta / Zipf) corollary. Confirms the efficiency collapse survives
#   discreteness -- the zeta h-index stays in the Gaussian regime -- and plots
#   the simulated variance ratio against the finite-n prediction and the
#   leading rate. The two deterministic curves are drawn on a dense grid to
#   n = 1e6 so the approach to the leading-rate asymptote is visible: the zeta
#   amplification carries O(1) terms (the log c shift and the scale derivative
#   c'(s)/c(s)) that vanish for the Pareto (c = 1), so its asymptote is
#   approached more slowly.
#
# Produces : figures/are_zeta.png; console output
# Reads    : nothing (self-contained simulation)
# Requires : VGAM, ggplot2  (CRAN packages)
# Run      : from the supplement root --  Rscript zeta_corollary_sim.R
# Seed     : set.seed(1729) (all stochastic scripts use this seed)
# ======================================================================

for (d_ in c("figures")) if (!dir.exists(d_)) dir.create(d_, showWarnings=FALSE)

# Zeta(s): P(X=k)=k^{-s}/zeta(s).  Tail index alpha=s-1 (regularly varying, Gaussian regime).
# Claim: Var(h_hat)/Var(H_n) ~ [c^{1/(1+alpha)}/((1+alpha)^2 V_s)] n^{-alpha/(1+alpha)} log^2 n -> 0,
#   V_s = Var_s(log X) = (log zeta)''(s).  Same rate as the continuous Pareto.
# VGAM: rzeta(n, shape=alpha) has exponent s=alpha+1.
suppressMessages(library(VGAM))
set.seed(1729)

# ---- zeta pieces via VGAM on a grid ----
sgrid <- seq(1.15, 8, by = 0.005)
zeta0 <- zeta(sgrid); zeta1 <- zeta(sgrid, deriv = 1); zeta2 <- zeta(sgrid, deriv = 2)
m_of_s <- -zeta1 / zeta0                       # E[log X] = -zeta'/zeta
s_of_m <- approxfun(m_of_s, sgrid, rule = 2)   # invert: s from mean(log X)
zeta_of_s <- approxfun(sgrid, zeta0, rule = 2)
V_of_s <- approxfun(sgrid, zeta2/zeta0 - (zeta1/zeta0)^2, rule = 2)  # Var(log X)

tail_em <- function(h, s) {                    # sum_{k>=h} k^-s, Euler-Maclaurin (h not nec. integer)
  h <- pmax(h, 1)
  h^(1 - s)/(s - 1) + 0.5 * h^(-s) + (s/12) * h^(-s - 1)
}
surv <- function(h, s) tail_em(h, s) / zeta_of_s(s)      # S_s(h)
pop_h <- function(n, s) uniroot(function(h) n*surv(h,s) - h, c(1, n))$root
hindex <- function(x){ xs<-sort(x,decreasing=TRUE); w<-which(xs>=seq_along(xs)); if(length(w)) max(w) else 0 }

run <- function(s, ns, B = 1500) {
  alpha <- s - 1; Vs <- V_of_s(s)
  # integer survival S^-(h)=P(X>=h) for a deterministic Var(H_n) (no Monte Carlo noise),
  # used in the finite-n prediction below.
  Kx <- 5e6; ksx <- (1:Kx)^(-s); Zx <- sum(ksx); csx <- c(0, cumsum(ksx))
  Sint <- function(h) (Zx - csx[h]) / Zx
  exactVarH <- function(n){ h <- 1:min(n, 6000); pr <- pbinom(h-1L, n, Sint(h), lower.tail=FALSE)
                            EH <- sum(pr); sum((2*h - 1)*pr) - EH^2 }
  cat(sprintf("\nZeta s=%.2f (alpha=%.2f)  V_s=Var(logX)=%.4f\n", s, alpha, Vs))
  cat(sprintf("%6s %8s %8s %10s %10s %9s %9s %9s\n","n","h_n","E[H_n]","Var(H_n)","Var(hhat)","ratio","rate(lead)","finite-n"))
  out <- data.frame()
  for (n in ns) {
    hn <- pop_h(n, s)
    cc <- 1/((s-1)*zeta_of_s(s))                  # tail constant c = 1/((s-1) zeta(s))
    H <- numeric(B); Hh <- numeric(B)
    for (b in 1:B) {
      x <- rzeta(n, shape = alpha)
      H[b] <- hindex(x)
      sh <- s_of_m(mean(log(x)))
      Hh[b] <- tryCatch(uniroot(function(h) n*surv(h,sh)-h, c(1, n))$root, error=function(e) NA)
    }
    vH <- var(H); vHh <- var(Hh, na.rm=TRUE)
    # asymptotic constant on the pure n-power includes the c^{1/(1+alpha)} factor
    ratio_th <- (cc^(1/(1+alpha))/((1+alpha)^2 * Vs)) * n^(-alpha/(1+alpha)) * log(n)^2
    # finite-n prediction: delta-method Var(hhat)=(dh/ds)^2/(n V_s), where dh/ds keeps the
    # log n + O(1) amplification (not just the leading log n), over the binomial-law Var(H_n).
    dh <- (pop_h(n, s + 1e-4) - pop_h(n, s - 1e-4)) / 2e-4
    ratio_ref <- (dh^2 / (n * Vs)) / exactVarH(n)
    cat(sprintf("%6d %8.2f %8.2f %10.3f %10.3f %9.4f %9.4f %9.4f\n",
                n, hn, mean(H), vH, vHh, vHh/vH, ratio_th, ratio_ref))
    out <- rbind(out, data.frame(n=n, hn=hn, EH=mean(H), VarH=vH, Varhh=vHh,
                                 ratio_sim=vHh/vH, ratio_th=ratio_th, ratio_ref=ratio_ref))
  }
  invisible(out)
}

res1 <- run(2.0, c(100,300,1000,3000), B=8000)   # alpha=1
res2 <- run(2.5, c(100,300,1000,3000), B=6000)   # alpha=1.5

# ---- deterministic curves on a dense n-grid to n=1e6 (no simulation involved) ----
s <- 2.0; alpha <- s - 1; Vs <- V_of_s(s); cc <- 1/((s-1)*zeta_of_s(s))
Kx <- 5e6; ksx <- (1:Kx)^(-s); Zx <- sum(ksx); csx <- c(0, cumsum(ksx))
Sint <- function(h) (Zx - csx[h]) / Zx
# pbinom evaluated only within 12 standard deviations of the mean; beyond that
# the probabilities are numerically 0/1 and far-tail evaluation is unstable.
p_ge <- function(k, n, p){
  mu <- n*p; sd <- sqrt(pmax(n*p*(1-p), 1e-12)); z <- (k-1-mu)/sd
  out <- numeric(length(k)); out[z < -12] <- 1
  mid <- abs(z) <= 12
  out[mid] <- pbinom(k[mid]-1L, n, p[mid], lower.tail=FALSE)
  out
}
exactVarH2 <- function(n){ h <- 1:min(n, 6000); pr <- p_ge(h, n, Sint(h))
                           EH <- sum(pr); sum((2*h - 1)*pr) - EH^2 }
ngrid <- unique(round(10^seq(2, 6, by = 0.1)))
grid <- data.frame(n = ngrid)
grid$ratio_th  <- (cc^(1/(1+alpha))/((1+alpha)^2 * Vs)) * grid$n^(-alpha/(1+alpha)) * log(grid$n)^2
grid$ratio_ref <- sapply(ngrid, function(n){
  dh <- (pop_h(n, s + 1e-4) - pop_h(n, s - 1e-4)) / 2e-4
  (dh^2 / (n * Vs)) / exactVarH2(n)
})

suppressMessages(library(ggplot2))
df <- res1
p <- ggplot() +
  geom_line(data=grid, aes(n, ratio_th, colour="Leading rate (Corollary 3)"), linewidth=0.8, linetype=2) +
  geom_line(data=grid, aes(n, ratio_ref, colour="Finite-n prediction"), linewidth=0.9) +
  geom_point(data=df, aes(n, ratio_sim, colour="Simulation"), size=2.2) +
  scale_x_log10(breaks=10^(2:6), labels=expression(10^2,10^3,10^4,10^5,10^6)) +
  scale_colour_manual(values=c("Leading rate (Corollary 3)"="grey55","Finite-n prediction"="#0072B2",
                               "Simulation"="#1A1A1A"), name=NULL,
                      breaks=c("Simulation","Finite-n prediction","Leading rate (Corollary 3)")) +
  guides(colour=guide_legend(override.aes=list(
    linetype=c(NA,1,2), shape=c(16,NA,NA), linewidth=c(NA,0.9,0.8)))) +
  labs(x="n (papers)", y=expression(Var(hat(h)[n])/Var(H[n]))) +
  theme_minimal(base_size=11) + theme(legend.position="top",
                                      plot.margin=margin(6, 8, 6, 8))
ggsave("figures/are_zeta.png", p, width=6.2, height=4.0, dpi=600)
cat("\nWrote figures/are_zeta.png\nDONE.\n")
