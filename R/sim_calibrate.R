# ======================================================================
# sim_calibrate.R
#
# Purpose  :
#   Solves the pre-specified calibration targets of the simulation truths.
#   Targets (medians over ASN records): h_n(F)/n_ref and the tail spread
#   q90/q50, where h_n is the root of n_ref S~(h) = h and q_p solves
#   S~(q) = 1 - p, S~ the log-linear extension of S^-.
#     short: n_ref = 16, h/n = 0.44, q90/q50 = 3.8
#     long : n_ref = 64, h/n = 0.20, q90/q50 = 3.8
#   Two-parameter families (C1 with c1 = 1.5, C2, T1, T2, T4) are solved by
#   nested root-finding on log-parameters: the scale-type parameter matches
#   h/n given the shape, the shape matches q90/q50. T3 takes the T1 shape
#   (sd) of the same setting, adds structural zeros with probability 0.3 and
#   re-solves mu for h/n only; T5 solves s for h/n only; T7 shares the T1
#   marginal. T2: base Burr(1.5, c2, s) with weight 0.9 plus a fixed heavy
#   component Burr(1.5, 2/3, 5 s) with weight 0.1 (tail index 1); q90/q50 is
#   unattainable with this component (minimum about 4.16), so the base c2 is
#   fixed at the C1 value of the same setting and s is solved for h/n only
#   (fallback = TRUE by construction). For the other two-parameter families an
#   unattainable q90/q50 would fall back to the shape closest to the target
#   (fallback = TRUE). Reports the implied
#   P(X = 0), q50, q90, q90/q50 and tail index; checks that every targeted
#   quantity is met to 1e-6 and writes the deterministic cell table.
#
# Produces : data/sim/calibration.csv, data/sim/cells.csv; console table
# Reads    : R/sim_common.R
# Requires : base R
# Run      : from the repository root --  Rscript R/sim_calibrate.R
# Seed     : none (deterministic)
# ======================================================================

source("R/sim_common.R")
for (d_ in c("data", "data/sim")) if (!dir.exists(d_)) dir.create(d_, showWarnings = FALSE)

TOL_ROOT <- 1e-11                     # absolute, on log-parameters (mu for PLN-type)

hn_ratio <- function(fam, p, n_ref) unname(truth_target(fam, p, n_ref)["h"]) / n_ref
qratio   <- function(fam, p) { q <- truth_quant(fam, p); q[2] / q[1] }

# Per family: shape range (log scale) and the law as a function of (shape, scale-type v)
cal_map <- list(
  C1 = list(range = log(c(0.05, 20)), vr = c(-15, 30),
            law = function(sh, v) list(c1 = C1_SHAPE1, c2 = sh, s = exp(v))),
  C2 = list(range = log(c(0.01, 200)), vr = c(-8, 20),
            law = function(sh, v) list(mu = exp(v), r = sh)),
  T1 = list(range = log(c(0.05, 5)), vr = c(-20, 30),
            law = function(sh, v) list(mu = v, sd = sh)),
  T2 = list(range = log(c(0.05, 50)), vr = c(-15, 30),
            law = function(sh, v) list(c1 = C1_SHAPE1, c2 = sh, s = exp(v))),
  T4 = list(range = log(c(0.05, 5)), vr = c(-20, 30),
            law = function(sh, v) list(mu = v, sd = sh)),
  T3 = list(vr = c(-20, 30), law = function(sh, v) list(mu = v, sd = sh)),
  T5 = list(vr = log(c(1e-3, 1e7)), law = function(sh, v) list(s = exp(v), alpha = T5_ALPHA))
)

# scale-type parameter v solving h_n/n_ref = rho for a given shape (h/n increasing in v);
# errors if the target is not bracketed on the family's range (shape infeasible)
solve_scale <- function(fam, sh, tg) {
  cm <- cal_map[[fam]]
  f <- function(v) hn_ratio(fam, cm$law(sh, v), tg$n_ref) - tg$rho
  f1 <- f(cm$vr[1]); f2 <- f(cm$vr[2])
  if (!(f1 < 0 && f2 > 0)) stop("h/n target not bracketed for shape ", sh)
  v <- uniroot(f, cm$vr, f.lower = f1, f.upper = f2, tol = TOL_ROOT)$root
  cm$law(sh, v)
}

calibrate_one <- function(fam, setting, calib_T1 = NULL, calib_C1 = NULL) {
  tg <- SIM_TARGETS[SIM_TARGETS$setting == setting, ]
  cm <- cal_map[[fam]]
  nsc <- NA_integer_; fallback <- FALSE; q_targeted <- TRUE
  if (fam == "T3") {
    p <- solve_scale("T3", calib_T1$sd[calib_T1$setting == setting], tg); q_targeted <- FALSE
  } else if (fam == "T5") {
    p <- solve_scale("T5", NA, tg); q_targeted <- FALSE
  } else if (fam == "T2") {                            # stipulated: base c2 = C1 c2, s for h/n only
    fallback <- TRUE; q_targeted <- FALSE
    p <- solve_scale("T2", calib_C1$c2[calib_C1$setting == setting], tg)
  } else {
    g <- function(lsh) qratio(fam, solve_scale(fam, exp(lsh), tg)) - tg$qratio
    grid <- seq(cm$range[1], cm$range[2], length.out = 31)
    gv <- vapply(grid, function(v) tryCatch(g(v), error = function(e) NA_real_), 0)
    ok <- !is.na(gv)
    sc <- which(ok[-1] & ok[-length(ok)] & sign(gv[-1]) != sign(gv[-length(gv)]))
    nsc <- length(sc)
    if (nsc) {
      lsh <- uniroot(g, grid[c(sc[1], sc[1] + 1)], tol = TOL_ROOT)$root
    } else {                                           # shape closest to the q90/q50 target
      fallback <- TRUE
      j <- which.min(abs(gv)); lo <- grid[max(1, j - 1)]; hi <- grid[min(length(grid), j + 1)]
      lsh <- optimize(function(v) abs(g(v)), c(lo, hi), tol = TOL_ROOT)$minimum
    }
    p <- solve_scale(fam, exp(lsh), tg)
  }
  hr <- hn_ratio(fam, p, tg$n_ref); q <- truth_quant(fam, p)
  v <- setNames(rep(NA_real_, 7), c("c1", "c2", "s", "mu", "r", "sd", "alpha")); v[names(p)] <- unlist(p)
  data.frame(family = fam, setting = setting, n_ref = tg$n_ref, rho_target = tg$rho,
             qratio_target = if (q_targeted) tg$qratio else NA_real_, t(v),
             rho_achieved = hr, qratio_achieved = q[2] / q[1], q50 = q[1], q90 = q[2],
             p0_implied = truth_p0(fam, p), tail_index = truth_tail_index(fam, p),
             err_rho = abs(hr - tg$rho),
             err_qratio = if (q_targeted) abs(q[2] / q[1] - tg$qratio) else NA_real_,
             n_sign_changes = nsc, fallback = fallback, stringsAsFactors = FALSE)
}

t0 <- proc.time()[3]
rows <- list()
for (fam in c("C1", "C2", "T1", "T2", "T4", "T5")) for (st in SIM$settings) {
  cC1 <- if (length(rows)) do.call(rbind, rows) else NULL
  rows[[length(rows) + 1]] <- calibrate_one(fam, st, calib_C1 = cC1[cC1$family == "C1", ])
  cat(sprintf("%-3s %-5s done (%.1fs)\n", fam, st, proc.time()[3] - t0))
}
calib <- do.call(rbind, rows)
cT1 <- calib[calib$family == "T1", ]
for (st in SIM$settings) calib <- rbind(calib, calibrate_one("T3", st, cT1))

# T7: same marginal as T1; achieved values recomputed through the T7 code path
t7 <- cT1; t7$family <- "T7"
for (i in seq_len(nrow(t7))) {
  p <- list(mu = t7$mu[i], sd = t7$sd[i]); q <- truth_quant("T7", p)
  t7$rho_achieved[i] <- hn_ratio("T7", p, t7$n_ref[i])
  t7$qratio_achieved[i] <- q[2] / q[1]; t7$q50[i] <- q[1]; t7$q90[i] <- q[2]
  t7$p0_implied[i] <- truth_p0("T7", p)
  t7$err_rho[i] <- abs(t7$rho_achieved[i] - t7$rho_target[i])
  t7$err_qratio[i] <- abs(t7$qratio_achieved[i] - t7$qratio_target[i])
}
calib <- rbind(calib, t7)
calib$pass <- calib$err_rho < 1e-6 & (is.na(calib$err_qratio) | calib$fallback | calib$err_qratio < 1e-6)
calib <- calib[order(match(calib$family, SIM$families), match(calib$setting, SIM$settings)), ]
rownames(calib) <- NULL

write.csv(calib, "data/sim/calibration.csv", row.names = FALSE)
print(calib, digits = 6)
if (any(calib$fallback)) cat("FLAG: fallback calibration, q90/q50 not targeted (see header):",
                             paste(calib$family[calib$fallback], calib$setting[calib$fallback]), "\n")
if (!all(calib$pass)) stop("calibration check failed: targets not met to 1e-6")
cat("calibration check passed: all targeted quantities met to 1e-6\n")

cells <- build_cells(calib)
write.csv(cells, "data/sim/cells.csv", row.names = FALSE)
cat(sprintf("cells.csv: %d cells (%d balanced, %d unbalanced)\n", nrow(cells),
            sum(!cells$unbalanced), sum(cells$unbalanced)))
