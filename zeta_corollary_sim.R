# ======================================================================
# zeta_corollary_sim.R
#
# Purpose  :
#   Discrete (zeta / Zipf) corollary. Confirms the efficiency collapse survives
#   discreteness -- the zeta h-index stays in the Gaussian regime -- and plots
#   the simulated against the predicted variance ratio.
#
# Produces : figures/are_zeta.png; console output
# Reads    : nothing (self-contained simulation)
# Requires : VGAM, ggplot2  (CRAN packages)
# Run      : from the supplement root --  Rscript zeta_corollary_sim.R
# Seed     : set.seed(1729) (all stochastic scripts use this seed)
# ======================================================================

## ensure output folders exist (some zip extractors drop empty dirs)
for (d_ in c("data","figures")) if (!dir.exists(d_)) dir.create(d_, showWarnings=FALSE)

# =====================================================================
# Zeta / discrete-Pareto corollary to Theorem 1 (FAST version).
# Zeta(s): P(X=k)=k^{-s}/zeta(s).  Tail index alpha=s-1 (RV tail, GAUSSIAN regime).
# Claim: Var(h_hat)/Var(H_n) ~ [c^{1/(1+alpha)}/((1+alpha)^2 V_s)] n^{-alpha/(1+alpha)} log^2 n -> 0,
#   V_s = Var_s(log X) = (log zeta)''(s).  Same RATE as continuous Pareto.
# Speed: precompute MLE inversion m(s)=-zeta'/zeta on a grid; closed-form
#   Euler-Maclaurin survival tail (no per-replicate series sums).
# VGAM: rzeta(n, shape=alpha) has exponent s=alpha+1.
# =====================================================================
suppressMessages(library(VGAM))
set.seed(1729)

# ---- exact-ish pieces via VGAM on a grid ----
sgrid <- seq(1.15, 8, by = 0.005)
zeta0 <- zeta(sgrid); zeta1 <- zeta(sgrid, deriv = 1); zeta2 <- zeta(sgrid, deriv = 2)
m_of_s <- -zeta1 / zeta0                       # E[log X] = -zeta'/zeta, increasing? decreasing in s
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
  # exact integer survival S^-(h)=P(X>=h) for a DETERMINISTIC Var(H_n) (no MC noise),
  # used in the refined finite-n prediction line below.
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
    # asymptotic constant on the pure n-power INCLUDES the c^{1/(1+alpha)} factor
    ratio_th <- (cc^(1/(1+alpha))/((1+alpha)^2 * Vs)) * n^(-alpha/(1+alpha)) * log(n)^2
    # refined finite-n prediction: delta-method Var(hhat)=(dh/ds)^2/(n V_s) [dh/ds keeps the
    # log n + O(1) amplification, not just leading log n] over the EXACT Var(H_n). Tracks the
    # simulation to MC noise; the gap to ratio_th is the dropped O(1)/binomial-factor terms.
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

suppressMessages(library(ggplot2))
figdir <- "figures"
df <- res1
p <- ggplot(df, aes(n)) +
  geom_line(aes(y=ratio_th, colour="Corollary 1 leading rate"), linewidth=0.9, linetype=2) +
  geom_line(aes(y=ratio_ref, colour="Finite-n prediction"), linewidth=1) +
  geom_point(aes(y=ratio_sim, colour="Simulation"), size=2.8) +
  scale_x_log10() +
  scale_colour_manual(values=c("Corollary 1 leading rate"="#8C8C8C","Finite-n prediction"="#2D5A3D",
                               "Simulation"="#C2453E"), name=NULL,
                      breaks=c("Simulation","Finite-n prediction","Corollary 1 leading rate")) +
  labs(x="n (papers)", y=expression(Var(hat(h)[n])/Var(H[n])),
       title="Zeta citations: variance ratio by sample size",
       subtitle="Simulation, finite-n prediction, and leading-rate curve") +
  theme_minimal(base_size=12) + theme(legend.position="top",
                                      plot.margin=margin(6, 8, 8, 8))
ggsave(file.path(figdir,"are_zeta.png"), p, width=6.2, height=4.2, dpi=150)
cat("\nWrote figure:", file.path(figdir,"are_zeta.png"), "\nDONE.\n")
