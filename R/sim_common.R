# ======================================================================
# sim_common.R
#
# Purpose  :
#   Function library for the simulation study of model-based h-index
#   estimators under misspecified citation laws. Defines
#   (i)   the discrete truths: S^-(k) = P(X >= k) at integers, samplers,
#         unit perturbations of a calibrated centre;
#   (ii)  the target: log-linear continuous extension of S^-, the root h_n
#         of n S~(h) = h, the integer target h_int and threshold decisions;
#   (iii) discretized working models X = floor(Y) (Pareto, lognormal,
#         Burr XII) fitted by interval-probability MLE, zero-truncated or
#         untruncated, and the NB2 likelihood on all counts;
#   (iv)  fit-validity rules, plug-in roots, complete pooling (lognormal,
#         Burr), partial pooling (Burr, MAP under a donor-estimated normal
#         prior on the log shapes), the Hill-Weissman plug-in and the
#         operational estimator (fallback to H when undefined);
#   (v)   the cell table and one L'Ecuyer-CMRG stream per (cell, field).
#
# Produces : nothing (sourced by sim_calibrate.R, sim_run.R, sim_aggregate.R)
# Reads    : nothing
# Requires : base R (stats, parallel); statmod optional (Gauss-Hermite
#            nodes, used only by the quadrature cross-check gh_nodes())
# Run      : source("R/sim_common.R") from the repository root
# Seed     : no draws at source time; per-(cell, field) streams from
#            sim_stream_seed() (L'Ecuyer-CMRG, master seed 20260930)
# ======================================================================

# ---------------------------------------------------------------------------
# Design constants
# ---------------------------------------------------------------------------
SIM <- list(
  master_seed   = 20260930L,
  U             = 30L,                 # units per field
  F_block       = 150L,                # fields per block
  F_max         = 600L,                # maximum fields per cell (4 blocks)
  n_grid        = c(12, 16, 25, 32, 64, 128, 256, 512, 1024),
  families      = c("C1", "C2", "T1", "T2", "T3", "T4", "T5", "T7"),
  settings      = c("short", "long"),
  taus          = c(0.1, 0.3),
  unbal_n       = c(16, 64, 256),
  unbal_tau     = 0.3,
  unbal_min_n   = 8,
  kappas        = c(2, 4),
  bound_margin  = 0.05,                # log-parameter distance to a box bound
  pp_min_donors = 10L,
  pp_eig_floor  = 1e-4,
  sens_untrunc  = TRUE,                # also fit untruncated LN and Burr
  mcse_target   = 0.02,
  ratio_band    = c(0.8, 1.25),
  quad_rel_tol  = 1e-12
)

# Calibration targets: h_n/n_ref and the tail spread q90/q50, with q_p solving
# S~(q) = 1 - p (log-linear extension). h/n = 0.44 and q90/q50 = 3.8 are the
# medians over ASN records with at least one work in the window, records built at
# the application deadlines (data/asn/derived/asn_records_summary.tsv).
SIM_TARGETS <- data.frame(setting = c("short", "long"), n_ref = c(16, 64),
                          rho = c(0.44, 0.20), qratio = c(3.8, 3.8),
                          stringsAsFactors = FALSE)

C1_SHAPE1 <- 1.5                     # Burr c1 of the C1 and T2 centres
T2_W      <- 0.1                     # weight of the heavy component
T2_HC1    <- 1.5                     # heavy component: fixed c1
T2_HC2    <- 2 / 3                   # heavy component: fixed c2 (tail index c1 c2 = 1)
T2_SF     <- 5                       # heavy component: scale 5 s (s = base scale)
T3_PI0    <- 0.30                    # structural-zero probability (T3 = T1 shape + zeros)
T4_AMAX   <- 15                      # paper age uniform on 1..15
T4_G      <- function(a) 1 - exp(-a / 4)
T5_ALPHA  <- 1.2
T5_BETA   <- 0.3
T7_LINE   <- 5                       # papers per research line
T7_SHARE  <- 0.3                     # between-line share of Var(log lambda)

BIG <- 1e100                         # objective value returned off-support

EST_FITS <- c("Pareto", "LN", "Burr", "NB2", "LNpool", "Burrpool", "BurrPP",
              if (SIM$sens_untrunc) c("LN_u", "Burr_u"))
EST_HW   <- paste0("HW", SIM$kappas)
EST_ALL  <- c("H", EST_FITS, EST_HW)
EST_PRIMARY <- setdiff(EST_ALL, c("LN_u", "Burr_u"))

# ---------------------------------------------------------------------------
# Numerics
# ---------------------------------------------------------------------------
log1pexp <- function(u) {                    # log(1 + exp(u))
  out <- u
  i <- which(!(u > 35))
  out[i] <- log1p(exp(u[i]))
  out
}

log1mexp <- function(d) {                    # log(1 - exp(-d)), d >= 0
  out <- rep(-Inf, length(d))
  out[is.na(d)] <- NaN
  i <- which(d > 0 & d <= log(2)); out[i] <- log(-expm1(-d[i]))
  j <- which(d > log(2));          out[j] <- log1p(-exp(-d[j]))
  out
}

lse2 <- function(a, b) {                     # log(exp(a) + exp(b))
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  out[m == -Inf] <- -Inf
  out
}

gh_nodes <- function(N) {                    # probabilists' Gauss-Hermite: E f(Z), Z ~ N(0,1)
  if (requireNamespace("statmod", quietly = TRUE)) {
    q <- statmod::gauss.quad.prob(N, dist = "normal")
    return(list(x = q$nodes, w = q$weights))
  }
  i <- seq_len(N - 1)
  J <- matrix(0, N, N)
  J[cbind(i, i + 1)] <- sqrt(i); J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}

# ---------------------------------------------------------------------------
# Continuous extension of S^- and the root of N S~(h) = h
#   logSm(k): log S^-(k) at integer k (vectorized or scalar), logSm(0) = 0.
#   S~(h) = S^-(k)^(k+1-h) S^-(k+1)^(h-k), k = floor(h).
#   Returns c(h = root, h_int = max{k : N S^-(k) >= k}); kmax is an integer
#   known (or expected) to satisfy N S^-(kmax + 1) < kmax + 1.
# ---------------------------------------------------------------------------
h_root <- function(logSm, logN, kmax) {
  ok <- function(k) isTRUE((logN + logSm(k)) >= log(k))
  lo <- 0; hi <- max(1, floor(kmax) + 1)
  while (ok(hi)) { lo <- hi; hi <- 2 * hi }
  while (hi - lo > 1) {
    mid <- (lo + hi) %/% 2
    if (ok(mid)) lo <- mid else hi <- mid
  }
  a <- if (lo == 0) 0 else logSm(lo)
  b <- logSm(lo + 1)
  if (!is.finite(a) || is.na(b)) return(c(h = NA_real_, h_int = lo))
  if (b == -Inf) return(c(h = lo, h_int = lo))
  phi <- function(h) logN + a + (h - lo) * (b - a) - log(h)
  low <- if (lo == 0) 1e-300 else lo
  r <- uniroot(phi, c(low, lo + 1), tol = 1e-12)$root
  c(h = r, h_int = lo)
}

# q_p solving S~(q) = 1 - p, 0 < p < 1 (continuous and piecewise smooth in the law)
q_tilde <- function(logSm, p) {
  L <- log1p(-p)
  ok <- function(k) isTRUE(logSm(k) >= L)
  lo <- 0; hi <- 1
  while (ok(hi)) { lo <- hi; hi <- 2 * hi }
  while (hi - lo > 1) {
    mid <- (lo + hi) %/% 2
    if (ok(mid)) lo <- mid else hi <- mid
  }
  a <- if (lo == 0) 0 else logSm(lo)
  b <- logSm(lo + 1)
  if (b == -Inf) return(lo)
  lo + (L - a) / (b - a)
}

S_tilde <- function(logSm, h) {              # continuous extension at real h >= 0
  k <- floor(h)
  exp((k + 1 - h) * logSm(k) + (h - k) * logSm(k + 1))
}

# ---------------------------------------------------------------------------
# Truths. Parameters per family (named list):
#   C1, T2: c1, c2, s     C2: mu (mean), r (size)
#   T1, T3, T4, T7: mu, sd (of log lambda / log Q)     T5: s, alpha
# ---------------------------------------------------------------------------
burr_logS <- function(y, c1, c2, s) {        # log P(Y >= y), Y ~ Burr XII
  out <- numeric(length(y))
  i <- y > 0
  out[i] <- -c2 * log1pexp(c1 * (log(y[i]) - log(s)))
  out
}

nb_logSm <- function(k, r, mu) {             # log P(X >= k), X ~ NB2(mu, r)
  out <- numeric(length(k))
  i <- k >= 1
  out[i] <- pnbinom(k[i] - 1, size = r, mu = mu, lower.tail = FALSE, log.p = TRUE)
  out
}

integ_safe <- function(f, a, b) {           # relative tolerance first; absolute floor if the segment is ~0
  tryCatch(integrate(f, a, b, rel.tol = SIM$quad_rel_tol, abs.tol = 0, subdivisions = 1000L)$value,
           error = function(e) integrate(f, a, b, rel.tol = 1e-10, abs.tol = 1e-300,
                                         subdivisions = 1000L, stop.on.error = FALSE)$value)
}

pln_S1 <- function(k, mu, sd) {              # P(X >= k), k >= 1 scalar; split at lambda = k and z = 0
  f <- function(z) dnorm(z) * ppois(k - 1, exp(mu + sd * z), lower.tail = FALSE)
  zs <- (log(k) - mu) / sd
  br <- sort(unique(c(-Inf, min(0, zs), max(0, zs), Inf)))
  tot <- 0
  for (j in seq_len(length(br) - 1)) tot <- tot + integ_safe(f, br[j], br[j + 1])
  max(tot, 0)
}

pln_p0 <- function(mu, sd) {                 # P(X = 0) = E exp(-lambda), split at lambda = 1
  f <- function(z) dnorm(z) * exp(-exp(mu + sd * z))
  zs <- -mu / sd
  br <- sort(unique(c(-Inf, min(0, zs), max(0, zs), Inf)))
  tot <- 0
  for (j in seq_len(length(br) - 1)) tot <- tot + integ_safe(f, br[j], br[j + 1])
  tot
}

pln_logSm <- function(k, mu, sd)
  vapply(k, function(kk) if (kk <= 0) 0 else log(pln_S1(kk, mu, sd)), 0)

truth_logSm <- function(fam, p, k) {
  switch(fam,
    C1 = burr_logS(k, p$c1, p$c2, p$s),
    T2 = {
      out <- lse2(log(1 - T2_W) + burr_logS(k, p$c1, p$c2, p$s),
                  log(T2_W) + burr_logS(k, T2_HC1, T2_HC2, T2_SF * p$s))
      out[k <= 0] <- 0
      out
    },
    C2 = nb_logSm(k, p$r, p$mu),
    T1 = , T7 = pln_logSm(k, p$mu, p$sd),
    T3 = {
      out <- log(1 - T3_PI0) + pln_logSm(k, p$mu, p$sd)
      out[k <= 0] <- 0
      out
    },
    T4 = vapply(k, function(kk) {
      if (kk <= 0) return(0)
      log(mean(vapply(seq_len(T4_AMAX), function(a)
        pln_S1(kk, p$mu + log(T4_G(a)), p$sd), 0)))
    }, 0),
    T5 = {
      l1 <- log1p(pmax(k, 0) / p$s)
      out <- log(0.5) + lse2(-p$alpha * l1, -(p$alpha + T5_BETA) * l1)
      out[k <= 0] <- 0
      out
    },
    stop("unknown family ", fam))
}

truth_p0 <- function(fam, p) {               # P(X = 0) = 1 - S^-(1)
  if (fam %in% c("T1", "T7")) return(pln_p0(p$mu, p$sd))
  if (fam == "T3") return(T3_PI0 + (1 - T3_PI0) * pln_p0(p$mu, p$sd))
  if (fam == "T4") return(mean(vapply(seq_len(T4_AMAX), function(a)
    pln_p0(p$mu + log(T4_G(a)), p$sd), 0)))
  if (fam == "C2") return(dnbinom(0, size = p$r, mu = p$mu))
  -expm1(truth_logSm(fam, p, 1))
}

rburr <- function(n, c1, c2, s) s * expm1(rexp(n) / c2)^(1 / c1)

truth_rx <- function(fam, p, n) {
  switch(fam,
    C1 = floor(rburr(n, p$c1, p$c2, p$s)),
    T2 = {
      hv <- runif(n) < T2_W
      floor(rburr(n, ifelse(hv, T2_HC1, p$c1), ifelse(hv, T2_HC2, p$c2), ifelse(hv, T2_SF * p$s, p$s)))
    },
    C2 = as.numeric(rnbinom(n, size = p$r, mu = p$mu)),
    T1 = as.numeric(rpois(n, exp(rnorm(n, p$mu, p$sd)))),
    T3 = {
      z <- runif(n) < T3_PI0
      x <- as.numeric(rpois(n, exp(rnorm(n, p$mu, p$sd))))
      x[z] <- 0
      x
    },
    T4 = {
      a <- sample.int(T4_AMAX, n, replace = TRUE)
      q <- exp(rnorm(n, p$mu, p$sd))
      as.numeric(rpois(n, q * T4_G(a)))
    },
    T5 = {
      hv <- runif(n) < 0.5
      a <- ifelse(hv, p$alpha + T5_BETA, p$alpha)
      floor(p$s * expm1(rexp(n) / a))
    },
    T7 = {
      nl <- ceiling(n / T7_LINE)
      u <- rnorm(nl, 0, sqrt(T7_SHARE) * p$sd)
      e <- rnorm(n, 0, sqrt(1 - T7_SHARE) * p$sd)
      line <- rep(seq_len(nl), each = T7_LINE)[seq_len(n)]
      as.numeric(rpois(n, exp(p$mu + u[line] + e)))
    },
    stop("unknown family ", fam))
}

# Unit law from the calibrated centre; z = (z1 scale, za first shape, zb second shape)
unit_law <- function(fam, cen, tau, z) {
  switch(fam,
    C1 = , T2 = list(c1 = cen$c1 * exp(tau * z[2]), c2 = cen$c2 * exp(tau * z[3]),
                     s = cen$s * exp(0.5 * z[1])),
    C2 = list(mu = cen$mu * exp(0.5 * z[1]), r = cen$r * exp(tau * z[2])),
    T1 = , T3 = , T4 = , T7 = list(mu = cen$mu + 0.5 * z[1], sd = cen$sd * exp(tau * z[2])),
    T5 = list(s = cen$s * exp(0.5 * z[1]), alpha = cen$alpha * exp(tau * z[2])),
    stop("unknown family ", fam))
}

truth_target <- function(fam, p, n)          # c(h, h_int) for a law at record size n
  h_root(function(k) truth_logSm(fam, p, k), log(n), n)

truth_quant <- function(fam, p, prob = c(0.5, 0.9))
  vapply(prob, function(pp) q_tilde(function(k) truth_logSm(fam, p, k), pp), 0)

truth_tail_index <- function(fam, p) {       # regular-variation index of S^- (NA if none)
  switch(fam, C1 = p$c1 * p$c2, T2 = min(T2_HC1 * T2_HC2, p$c1 * p$c2), T5 = p$alpha, NA_real_)
}

# ---------------------------------------------------------------------------
# Working models on the optimization scale th
#   Pareto: log alpha            (S(y) = y^-alpha, y >= 1)
#   LN    : mu, log sigma        (S(y) = Phi((mu - log y)/sigma))
#   Burr  : log c1, log c2, log s (S(y) = {1 + (y/s)^c1}^-c2)
#   NB2   : log r, log mu
# ---------------------------------------------------------------------------
WM_BOUNDS <- list(
  Pareto = list(lower = log(0.05), upper = log(20)),
  LN     = list(lower = c(-10, log(0.05)), upper = c(15, log(10))),
  Burr   = list(lower = log(c(0.05, 0.01, 1e-3)), upper = log(c(20, 50, 1e6))),
  NB2    = list(lower = log(c(0.01, 1e-3)), upper = log(c(1e3, 1e5)))
)

wm_logS <- function(fam, y, th, grad = FALSE) {   # log S(y) and d/dth at y >= 0
  ny <- length(y); pos <- y > 0
  if (fam == "Pareto") {
    ls <- -exp(th[1]) * log(pmax(y, 1))
    return(list(ls = ls, g = if (grad) matrix(ls, ny, 1)))
  }
  if (fam == "LN") {
    ls <- numeric(ny); g <- if (grad) matrix(0, ny, 2)
    if (any(pos)) {
      sg <- exp(th[2]); z <- (th[1] - log(y[pos])) / sg
      lz <- pnorm(z, log.p = TRUE); ls[pos] <- lz
      if (grad) { rt <- exp(dnorm(z, log = TRUE) - lz); g[pos, 1] <- rt / sg; g[pos, 2] <- -rt * z }
    }
    return(list(ls = ls, g = g))
  }
  if (fam == "Burr") {
    ls <- numeric(ny); g <- if (grad) matrix(0, ny, 3)
    if (any(pos)) {
      c1 <- exp(th[1]); c2 <- exp(th[2])
      u <- c1 * (log(y[pos]) - th[3]); L <- log1pexp(u)
      ls[pos] <- -c2 * L
      if (grad) { pu <- plogis(u); g[pos, 1] <- -c2 * pu * u; g[pos, 2] <- -c2 * L; g[pos, 3] <- c2 * pu * c1 }
    }
    return(list(ls = ls, g = g))
  }
  stop("wm_logS: unknown family ", fam)
}

wm_logSm <- function(fam, k, th) {           # fitted log S^-(k) at integers
  if (fam == "NB2") return(nb_logSm(k, exp(th[1]), exp(th[2])))
  wm_logS(fam, k, th)$ls
}

wm_prep <- function(x, trunc) {
  xx <- if (trunc) x[x >= 1] else x
  kv <- sort(unique(xx))
  list(kv = kv, wv = tabulate(match(xx, kv), length(kv)), m = length(xx), n = length(x),
       xx = xx)
}

GL4_X <- c(-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526)
GL4_W <- c(0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538) / 2

mills <- function(z) exp(dnorm(z, log = TRUE) - pnorm(z, log.p = TRUE))   # phi(z)/Phi(z)

# log P(k <= Y < k+1) = log S(k) + log(1 - exp(-d)), d = log S(k) - log S(k+1),
# with d computed without subtracting the two logs (no cancellation at large k),
# and its gradient in th: G = dlogS(k) + dd / expm1(d).
wm_interval <- function(fam, k, th) {
  A <- wm_logS(fam, k, th, grad = TRUE)
  nk <- length(k); np <- length(th)
  d <- numeric(nk); dd <- matrix(0, nk, np)
  z0 <- which(k == 0); kp <- which(k > 0)
  if (length(z0)) {                                   # S(0) = 1: d = -log S(1)
    B1 <- wm_logS(fam, 1, th, grad = TRUE)
    d[z0] <- -B1$ls; dd[z0, ] <- matrix(-B1$g[1, ], length(z0), np, byrow = TRUE)
  }
  if (length(kp)) {
    kk <- k[kp]; l <- log1p(1 / kk)
    if (fam == "Pareto") {
      d[kp] <- exp(th[1]) * l; dd[kp, 1] <- d[kp]
    } else if (fam == "Burr") {
      c1 <- exp(th[1]); c2 <- exp(th[2])
      u <- c1 * (log(kk) - th[3]); p <- plogis(u); pq <- p * plogis(-u); E <- expm1(c1 * l)
      q <- p * E; d[kp] <- c2 * log1p(q)
      dd[kp, 1] <- c2 * (pq * u * E + p * c1 * l * (E + 1)) / (1 + q)
      dd[kp, 2] <- d[kp]
      dd[kp, 3] <- c2 * (-c1 * pq * E) / (1 + q)
    } else if (fam == "LN") {
      sg <- exp(th[2]); zk <- (th[1] - log(kk)) / sg; w <- l / sg; zk1 <- zk - w
      nar <- w < 0.05
      dv <- d1 <- d2 <- numeric(length(kk))
      if (any(!nar)) {                                # wide interval: direct difference
        i <- which(!nar); la <- mills(zk[i]); lb <- mills(zk1[i])
        dv[i] <- pnorm(zk[i], log.p = TRUE) - pnorm(zk1[i], log.p = TRUE)
        d1[i] <- (la - lb) / sg; d2[i] <- -la * zk[i] + lb * zk1[i]
      }
      if (any(nar)) {                                 # narrow interval: d = int_{zk1}^{zk} mills(t) dt, 4-pt Gauss-Legendre
        i <- which(nar); s0 <- s1 <- s2 <- 0
        for (j in 1:4) {
          t <- zk1[i] + w[i] * (1 + GL4_X[j]) / 2; lm <- mills(t); lmp <- -lm * (t + lm)
          s0 <- s0 + GL4_W[j] * lm; s1 <- s1 + GL4_W[j] * lmp; s2 <- s2 + GL4_W[j] * lmp * t
        }
        dv[i] <- w[i] * s0; d1[i] <- w[i] * s1 / sg; d2[i] <- -dv[i] - w[i] * s2
      }
      d[kp] <- dv; dd[kp, 1] <- d1; dd[kp, 2] <- d2
    } else stop("wm_interval: unknown family ", fam)
  }
  list(lp = A$ls + log1mexp(d), G = A$g + dd / expm1(d))
}

ic_nll <- function(fam, th, dat, trunc) {    # interval-probability negative log-likelihood
  I <- wm_interval(fam, dat$kv, th)
  if (any(!is.finite(I$lp))) return(list(f = Inf, g = rep(NA_real_, length(th))))
  f <- -sum(dat$wv * I$lp); g <- -colSums(dat$wv * I$G)
  if (trunc) { S1 <- wm_logS(fam, 1, th, grad = TRUE); f <- f + dat$m * S1$ls; g <- g + dat$m * S1$g[1, ] }
  list(f = f, g = g)
}

nb_nll <- function(th, dat) {                # NB2 negative log-likelihood on all counts
  r <- exp(th[1]); mu <- exp(th[2]); k <- dat$kv; w <- dat$wv
  lp <- dnbinom(k, size = r, mu = mu, log = TRUE)
  dr <- digamma(k + r) - digamma(r) - log1p(mu / r) + (mu - k) / (r + mu)
  dm <- k / mu - (r + k) / (r + mu)
  list(f = -sum(w * lp), g = -c(sum(w * dr) * r, sum(w * dm) * mu))
}

wm_objective <- function(fam, dat, trunc, th_fixed, free, prior = NULL) {
  cache <- new.env()
  ev <- function(p) {
    if (!is.null(cache$p) && identical(p, cache$p)) return(cache$res)
    th <- th_fixed; th[free] <- p
    res <- if (fam == "NB2") nb_nll(th, dat) else ic_nll(fam, th, dat, trunc)
    if (!is.null(prior) && is.finite(res$f)) {
      dv <- th[prior$idx] - prior$mu; Sd <- drop(prior$Sinv %*% dv)
      res$f <- res$f + 0.5 * sum(dv * Sd)
      res$g[prior$idx] <- res$g[prior$idx] + Sd
    }
    if (!is.finite(res$f) || any(!is.finite(res$g))) res <- list(f = BIG, g = rep(0, length(th)))
    res$g <- res$g[free]
    cache$p <- p; cache$res <- res
    res
  }
  list(fn = function(p) ev(p)$f, gr = function(p) ev(p)$g)
}

wm_starts <- function(fam, dat) {            # >= 4 starts incl. a moment start (full th)
  xx <- dat$xx
  ly <- log(xx + 0.5)
  s_ly <- if (length(ly) > 1) sd(ly) else 0
  s_ly <- max(s_ly, 0.2)
  switch(fam,
    Pareto = {
      a0 <- 1 / max(mean(log(xx + 1)), 0.05)
      as.list(log(c(a0, 0.5, 1, 2)))
    },
    LN = list(c(mean(ly), log(s_ly)), c(mean(ly), log(2 * s_ly)),
              c(mean(ly) - s_ly, log(1.5 * s_ly)), c(median(ly), log(0.7 * s_ly))),
    Burr = {
      c10 <- pi / (sqrt(3) * s_ly); s0 <- median(ly)   # log-logistic (c2 = 1) moment start
      list(c(log(c10), 0, s0), c(log(2 * c10), log(0.5), s0),
           c(log(c10), log(2), s0 + log(2)), c(log(3), log(0.3), s0))
    },
    NB2 = {
      mu0 <- max(mean(xx), 0.01)
      v <- if (length(xx) > 1) var(xx) else 2 * mu0
      r0 <- min(max(mu0^2 / max(v - mu0, 0.05 * mu0), 0.02), 900)
      list(c(log(r0), log(mu0)), c(log(0.3 * r0), log(mu0)),
           c(log(3 * r0), log(mu0)), c(0, log(mu0)))
    })
}

# Bounded multi-start L-BFGS-B. fixed: full-length vector, NA = free.
# prior: list(idx, mu, Sinv) adds 0.5 (th[idx]-mu)' Sinv (th[idx]-mu).
# Validity: converged, finite objective, every free parameter >= bound_margin
# from its box bound. Exception (zero-truncated Burr): the scale at its LOWER
# bound is the Pareto limit of the conditional law (S(y)/S(1) -> y^(-c1 c2));
# such a fit is defined and flagged at_pl = TRUE; its Hessian check uses the
# block of the remaining free parameters.
fit_wm <- function(x, fam, trunc, fixed = NULL, prior = NULL, starts = NULL,
                   hess = FALSE) {
  B <- WM_BOUNDS[[fam]]; np <- length(B$lower)
  if (is.null(fixed)) fixed <- rep(NA_real_, np)
  free <- which(is.na(fixed))
  out <- list(th = rep(NA_real_, np), def = FALSE, conv = NA_integer_, value = NA_real_,
              bnd = NA, at_pl = FALSE, nok = 0L, pd = NA, hess = NULL)
  dat <- wm_prep(x, trunc)
  if (dat$m < 1) return(out)
  if (is.null(starts)) starts <- wm_starts(fam, dat)
  obj <- wm_objective(fam, dat, trunc, fixed, free, prior)
  lo <- B$lower[free]; up <- B$upper[free]
  best <- NULL
  for (st in starts) {
    p0 <- pmin(pmax(st[free], lo + 0.01), up - 0.01)
    o <- tryCatch(optim(p0, obj$fn, obj$gr, method = "L-BFGS-B", lower = lo, upper = up,
                        control = list(maxit = 1000)), error = function(e) NULL)
    if (is.null(o)) next
    if (o$convergence == 0 && o$value < BIG / 10) out$nok <- out$nok + 1L
    if (is.null(best) || o$value < best$value) best <- o
  }
  if (is.null(best)) return(out)
  if (best$convergence != 0 && best$value < BIG / 10) {   # one restart from the best point
    o <- tryCatch(optim(best$par, obj$fn, obj$gr, method = "L-BFGS-B", lower = lo, upper = up,
                        control = list(maxit = 1000)), error = function(e) NULL)
    if (!is.null(o) && o$value <= best$value + 1e-8 * (1 + abs(best$value))) best <- o
  }
  th <- fixed; th[free] <- best$par
  out$th <- th; out$conv <- best$convergence; out$value <- best$value
  low_hit <- best$par - lo < SIM$bound_margin
  up_hit  <- up - best$par < SIM$bound_margin
  out$bnd <- any(low_hit | up_hit)
  s_pos <- if (fam == "Burr" && trunc) match(3L, free) else NA_integer_
  out$at_pl <- !is.na(s_pos) && low_hit[s_pos]
  other_hit <- low_hit | up_hit
  if (out$at_pl) other_hit[s_pos] <- up_hit[s_pos]
  out$def <- best$convergence == 0 && best$value < BIG / 10 && !any(other_hit)
  if (hess && out$def) {
    H <- tryCatch(optimHess(best$par, obj$fn, obj$gr), error = function(e) NULL)
    keep <- if (out$at_pl) -s_pos else seq_along(free)
    Hk <- if (is.null(H)) NULL else H[keep, keep, drop = FALSE]
    out$pd <- !is.null(Hk) && all(is.finite(Hk)) && length(Hk) > 0 &&
      min(eigen(Hk, symmetric = TRUE, only.values = TRUE)$values) > 0
    out$hess <- H
  }
  out
}

# Plug-in root in the target's convention.
#   truncated  : m S~_th(h) / S_th(1) = h   (m = #{X >= 1})
#   untruncated: n S~_th(h) = h
plugin_wm <- function(fam, th, n, m, trunc) {
  lS <- function(k) wm_logSm(fam, k, th)
  r <- if (trunc) h_root(lS, log(m) - lS(1), m) else h_root(lS, log(n), n)
  unname(r["h"])
}

# ---------------------------------------------------------------------------
# Empirical h-index and Hill-Weissman plug-in
# ---------------------------------------------------------------------------
h_emp <- function(x) { xs <- sort(x, decreasing = TRUE); sum(xs >= seq_along(xs)) }

hw_est <- function(x, H, kappa) {
  k <- floor(kappa * H); m <- sum(x >= 1)
  if (k < 1 || k + 1 > m) return(NA_real_)
  xs <- sort(x, decreasing = TRUE); Xk1 <- xs[k + 1]
  if (Xk1 < 1) return(NA_real_)
  g <- mean(log(xs[seq_len(k)])) - log(Xk1)
  if (!(g > 0)) return(NA_real_)
  a <- 1 / g
  exp((log(k) + a * log(Xk1)) / (1 + a))
}

# ---------------------------------------------------------------------------
# Partial-pooling prior from donors: defined Burr fits with PD Hessian, excluding
# fits at the Pareto limit (a boundary MLE; its Hessian is not the curvature at a
# stationary point, so its inverse is not used as a sampling variance)
# ---------------------------------------------------------------------------
pp_prior <- function(dfits) {
  psi <- t(vapply(dfits, function(f) f$th[1:2], c(0, 0)))
  V <- lapply(dfits, function(f) solve(f$hess)[1:2, 1:2])
  mu <- colMeans(psi)
  S <- cov(psi) - Reduce(`+`, V) / length(V)
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  ev <- pmax(e$values, SIM$pp_eig_floor)
  list(idx = 1:2, mu = mu,
       Sigma = e$vectors %*% diag(ev) %*% t(e$vectors),
       Sinv = e$vectors %*% diag(1 / ev) %*% t(e$vectors))
}

# ---------------------------------------------------------------------------
# One field of U units (the RNG must already be positioned on the field's stream)
#   cell: one row of the cell table (family, setting, tau, n, unbalanced, centre)
# ---------------------------------------------------------------------------
cell_centre <- function(cell) {
  nm <- c("c1", "c2", "s", "mu", "r", "sd", "alpha")
  v <- as.list(unlist(cell[nm])); v[!is.na(unlist(v))]
}

sim_field <- function(cell, U = SIM$U) {
  fam <- as.character(cell$family); tau <- cell$tau; n <- cell$n
  cen <- cell_centre(cell)
  Z <- matrix(rnorm(4 * U), U, 4, byrow = TRUE)        # z1 scale, za, zb shapes, zn record size
  laws <- lapply(seq_len(U), function(u) unit_law(fam, cen, tau, Z[u, ]))
  n_u <- if (isTRUE(as.logical(cell$unbalanced)))
    pmax(SIM$unbal_min_n, round(n * exp(0.5 * Z[, 4]))) else rep(n, U)
  X <- lapply(seq_len(U), function(u) truth_rx(fam, laws[[u]], n_u[u]))

  tg <- t(vapply(seq_len(U), function(u) truth_target(fam, laws[[u]], n_u[u]), c(h = 0, h_int = 0)))
  m_u <- vapply(X, function(x) sum(x >= 1), 0)
  H <- vapply(X, h_emp, 0)

  fits <- lapply(seq_len(U), function(u) {
    x <- X[[u]]
    f <- list(Pareto = fit_wm(x, "Pareto", TRUE),
              LN     = fit_wm(x, "LN", TRUE),
              Burr   = fit_wm(x, "Burr", TRUE, hess = TRUE),
              NB2    = fit_wm(x, "NB2", FALSE))
    if (SIM$sens_untrunc) {
      f$LN_u   <- fit_wm(x, "LN", FALSE)
      f$Burr_u <- fit_wm(x, "Burr", FALSE)
    }
    f
  })
  defLN <- vapply(fits, function(f) f$LN$def, TRUE)
  defB  <- vapply(fits, function(f) f$Burr$def, TRUE)
  pdB   <- vapply(fits, function(f) isTRUE(f$Burr$pd), TRUE)
  plB   <- vapply(fits, function(f) isTRUE(f$Burr$at_pl), TRUE)

  ndon <- matrix(0L, U, 3, dimnames = list(NULL, c("LNpool", "Burrpool", "BurrPP")))
  for (u in seq_len(U)) {
    x <- X[[u]]; don <- setdiff(seq_len(U), u)
    dl <- don[defLN[don]]; db <- don[defB[don]]; dp <- don[defB[don] & pdB[don] & !plB[don]]
    ndon[u, ] <- c(length(dl), length(db), length(dp))
    na_fit <- list(th = rep(NA_real_, 3), def = FALSE, conv = NA_integer_, value = NA_real_,
                   bnd = NA, at_pl = FALSE, nok = 0L, pd = NA)
    # complete pooling: donor-median log shapes, scale fitted on the unit
    fits[[u]]$LNpool <- if (length(dl) >= 1) {
      lsg <- median(vapply(fits[dl], function(f) f$LN$th[2], 0))
      mu0 <- mean(log(x[x >= 1] + 0.5))
      fit_wm(x, "LN", TRUE, fixed = c(NA, lsg),
             starts = lapply(mu0 + c(0, -1, 1, 2), function(a) c(a, lsg)))
    } else na_fit
    fits[[u]]$Burrpool <- if (length(db) >= 1) {
      psi <- apply(vapply(fits[db], function(f) f$Burr$th[1:2], c(0, 0)), 1, median)
      s0 <- median(log(x[x >= 1] + 0.5))
      fit_wm(x, "Burr", TRUE, fixed = c(psi, NA),
             starts = lapply(s0 + c(0, -1, 1, 2), function(a) c(psi, a)))
    } else na_fit
    # partial pooling: MAP with normal prior on (log c1, log c2)
    fits[[u]]$BurrPP <- if (length(dp) >= SIM$pp_min_donors) {
      pr <- pp_prior(lapply(fits[dp], `[[`, "Burr"))
      s0 <- median(log(x[x >= 1] + 0.5))
      st <- list(c(pr$mu, s0), c(pr$mu, s0 + 1), c(pr$mu, s0 - 1))
      if (all(is.finite(fits[[u]]$Burr$th))) st <- c(list(fits[[u]]$Burr$th), st)
      if (all(is.finite(fits[[u]]$Burrpool$th))) st <- c(st, list(c(pr$mu, fits[[u]]$Burrpool$th[3])))
      if (length(st) < 4) st <- c(st, list(c(pr$mu, s0 + 2)))
      f <- fit_wm(x, "Burr", TRUE, prior = pr, starts = st, hess = TRUE)
      f$def <- f$def && isTRUE(f$pd)
      f
    } else na_fit
  }

  wm_fam <- c(Pareto = "Pareto", LN = "LN", Burr = "Burr", NB2 = "NB2", LNpool = "LN",
              Burrpool = "Burr", BurrPP = "Burr", LN_u = "LN", Burr_u = "Burr")
  wm_tr  <- c(Pareto = TRUE, LN = TRUE, Burr = TRUE, NB2 = FALSE, LNpool = TRUE,
              Burrpool = TRUE, BurrPP = TRUE, LN_u = FALSE, Burr_u = FALSE)

  est <- matrix(NA_real_, U, length(EST_ALL), dimnames = list(NULL, EST_ALL))
  est[, "H"] <- H
  for (u in seq_len(U)) {
    for (e in EST_FITS) {
      f <- fits[[u]][[e]]
      if (isTRUE(f$def)) {
        v <- tryCatch(plugin_wm(wm_fam[[e]], f$th, n_u[u], m_u[u], wm_tr[[e]]),
                      error = function(err) NA_real_)
        if (is.finite(v)) est[u, e] <- v
      }
    }
    for (j in seq_along(SIM$kappas)) est[u, EST_HW[j]] <- hw_est(X[[u]], H[u], SIM$kappas[j])
  }
  def <- !is.na(est)
  op <- est; op[!def] <- matrix(H, U, length(EST_ALL))[!def]

  # unit rows
  lawm <- t(vapply(laws, function(p) {
    v <- setNames(rep(NA_real_, 7), c("c1", "c2", "s", "mu", "r", "sd", "alpha"))
    v[names(p)] <- unlist(p); v }, setNames(numeric(7), c("c1", "c2", "s", "mu", "r", "sd", "alpha"))))
  colnames(lawm) <- paste0("law_", colnames(lawm))
  units <- data.frame(unit = seq_len(U), lawm, z1 = Z[, 1], za = Z[, 2], zb = Z[, 3], zn = Z[, 4],
                      n_u = n_u, m_u = m_u, h_n = tg[, "h"], h_int = tg[, "h_int"],
                      check.names = FALSE)
  units <- cbind(units,
                 setNames(as.data.frame(est), paste0("est_", EST_ALL)),
                 setNames(as.data.frame(def), paste0("def_", EST_ALL)),
                 setNames(as.data.frame(op),  paste0("op_",  EST_ALL)))
  for (e in EST_FITS) {
    th <- t(vapply(fits, function(f) { v <- f[[e]]$th; c(v, rep(NA, 3 - length(v))) }, numeric(3)))
    units[[paste0(e, "_th1")]] <- th[, 1]; units[[paste0(e, "_th2")]] <- th[, 2]
    units[[paste0(e, "_th3")]] <- th[, 3]
    units[[paste0(e, "_conv")]] <- vapply(fits, function(f) as.integer(f[[e]]$conv), 0L)
    units[[paste0(e, "_nll")]]  <- vapply(fits, function(f) as.numeric(f[[e]]$value), 0)
    units[[paste0(e, "_bnd")]]  <- vapply(fits, function(f) as.logical(f[[e]]$bnd), TRUE)
    units[[paste0(e, "_nok")]]  <- vapply(fits, function(f) as.integer(f[[e]]$nok), 0L)
    units[[paste0(e, "_fitdef")]] <- vapply(fits, function(f) isTRUE(f[[e]]$def), TRUE)
  }
  units$Burr_pd   <- pdB
  units$Burr_at_pareto_limit     <- plB
  units$Burrpool_at_pareto_limit <- vapply(fits, function(f) isTRUE(f$Burrpool$at_pl), TRUE)
  units$BurrPP_at_pareto_limit   <- vapply(fits, function(f) isTRUE(f$BurrPP$at_pl), TRUE)
  units$BurrPP_pd <- vapply(fits, function(f) as.logical(f$BurrPP$pd), TRUE)
  units$LNpool_ndon <- ndon[, "LNpool"]; units$Burrpool_ndon <- ndon[, "Burrpool"]
  units$BurrPP_ndon <- ndon[, "BurrPP"]

  # decision rows: t in {h_int - 2, ..., h_int + 2}, t >= 1
  dec <- do.call(rbind, lapply(seq_len(U), function(u) {
    tt <- (tg[u, "h_int"] - 2):(tg[u, "h_int"] + 2); tt <- tt[tt >= 1]
    if (!length(tt)) return(NULL)
    d <- data.frame(unit = u, t = tt, delta = tt <= tg[u, "h_int"], dist = tt - tg[u, "h"])
    for (e in EST_ALL) d[[paste0("dec_", e)]] <- op[u, e] >= tt
    d
  }))
  list(units = units, dec = dec)
}

# ---------------------------------------------------------------------------
# Cell table (deterministic order) and RNG streams
# ---------------------------------------------------------------------------
build_cells <- function(calib) {
  main <- expand.grid(n = SIM$n_grid, tau = SIM$taus, setting = SIM$settings,
                      family = SIM$families, stringsAsFactors = FALSE)
  main <- main[, c("family", "setting", "tau", "n")]
  main$unbalanced <- FALSE
  ub <- expand.grid(n = SIM$unbal_n, setting = SIM$settings, family = SIM$families,
                    stringsAsFactors = FALSE)
  ub <- data.frame(family = ub$family, setting = ub$setting, tau = SIM$unbal_tau, n = ub$n,
                   unbalanced = TRUE, stringsAsFactors = FALSE)
  cells <- rbind(main, ub)
  cells$cell <- seq_len(nrow(cells))
  key <- c("family", "setting")
  pc <- c("c1", "c2", "s", "mu", "r", "sd", "alpha")
  cells <- merge(cells, calib[, c(key, pc)], by = key, all.x = TRUE, sort = FALSE)
  cells <- cells[order(cells$cell), c("cell", "family", "setting", "tau", "n", "unbalanced", pc)]
  rownames(cells) <- NULL
  cells
}

# Code/design fingerprint stored in every raw output and checked on resume and
# aggregation: md5 of sim_common.R, sim_run.R and data/sim/cells.csv.
SIM_FP_FILES <- c("R/sim_common.R", "R/sim_run.R", "data/sim/cells.csv")
sim_fingerprint <- function() {
  if (!all(file.exists(SIM_FP_FILES))) stop("fingerprint files missing: ",
    paste(SIM_FP_FILES[!file.exists(SIM_FP_FILES)], collapse = ", "))
  paste(basename(SIM_FP_FILES), unname(tools::md5sum(SIM_FP_FILES)), sep = "=", collapse = ";")
}

# Stream J = (cell - 1) * F_max + field, obtained by J applications of
# parallel::nextRNGStream to the L'Ecuyer-CMRG state after set.seed(master_seed).
sim_stream_seed <- function(cell, field) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(SIM$master_seed)
  s <- .Random.seed
  J <- (cell - 1) * SIM$F_max + field
  for (i in seq_len(J)) s <- parallel::nextRNGStream(s)
  s
}
