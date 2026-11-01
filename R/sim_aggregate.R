# ======================================================================
# sim_aggregate.R
#
# Purpose  :
#   Aggregates the per-(cell, block) simulation output. Per cell and
#   estimator: MSE ratio against H for the operational estimator
#   (unconditional; primary) and on defined units (conditional;
#   diagnostic), with field-level delta-method MCSE
#   Z_b = (A_b - r B_b)/mean(B), MCSE = sd(Z)/sqrt(F), where A_b and B_b
#   are the field sums of squared errors of T and of H; the same ratio
#   against the integer target h_int (sensitivity); bias; definedness.
#   Stopping rule: a cell needs another block of 150 fields while any
#   estimator with ratio in [0.8, 1.25] has MCSE > 0.02 and F < 600.
#   Crossing n: per (family, setting, tau, estimator) on the balanced
#   grid, the n where the ratio crosses 1 (log-log interpolation).
#   Threshold decisions: P{1(T >= t) != delta_t} by cell and estimator,
#   overall, for |t - h_n| <= 1 and by |t - h_n| bin.
#
# Produces : <out_dir>/results.rds, <out_dir>/summary.csv,
#            <out_dir>/tasklist_next.txt (cell block pairs still needed),
#            <fig> (MSE ratio vs n, one panel per truth x setting x tau)
# Reads    : R/sim_common.R, <raw_dir>/cell_<i>_block_<b>.rds
# Requires : base R; ggplot2 for the figure (skipped if absent)
# Checks   : every raw file must carry the current fingerprint (md5 of
#            sim_common.R, sim_run.R, data/sim/cells.csv), else error. Each
#            cell must hold contiguous field ids 1..F in complete blocks of
#            150 with U units per field, and no estimator with ratio in
#            [0.8, 1.25] may have an NA MCSE; otherwise the script errors
#            before writing tasklist_next, unless the fourth argument is
#            "partial" (smoke/debug runs: the checks warn and no tasklist is
#            written).
# Run      : from the repository root --
#              Rscript R/sim_aggregate.R [raw_dir] [out_dir] [fig] [partial]
#            defaults: data/sim/raw  data/sim  figures/sim_mse_ratio_diagnostic.png (diagnostic; the paper figure is drawn by R/sim_report.R)
# Seed     : none (deterministic)
# ======================================================================

args <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args) >= 1) args[1] else "data/sim/raw"
out_dir <- if (length(args) >= 2) args[2] else "data/sim"
fig     <- if (length(args) >= 3) args[3] else "figures/sim_mse_ratio_diagnostic.png"
partial <- length(args) >= 4 && identical(args[4], "partial")
source("R/sim_common.R")
for (d_ in c(out_dir, dirname(fig))) if (!dir.exists(d_)) dir.create(d_, recursive = TRUE, showWarnings = FALSE)

files <- list.files(raw_dir, pattern = "^cell_[0-9]+_block_[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no raw files in ", raw_dir)
L <- lapply(files, readRDS)
fp <- sim_fingerprint()
bad <- vapply(L, function(o) !identical(o$meta$fingerprint, fp), TRUE)
if (any(bad)) stop(sprintf("%d raw file(s) carry a different fingerprint than the current code/cells (e.g. %s)\n  current: %s",
                           sum(bad), basename(files[which(bad)[1]]), fp))
units <- do.call(rbind, lapply(L, `[[`, "units"))
dec   <- do.call(rbind, lapply(L, `[[`, "dec"))
if (any(duplicated(units[, c("cell", "field", "unit")]))) stop("duplicated (cell, field, unit) rows")
ests <- intersect(EST_ALL, sub("^op_", "", grep("^op_", names(units), value = TRUE)))
keys <- c("cell", "family", "setting", "tau", "n", "unbalanced")

ratio_mcse <- function(A, B) {
  if (!length(A) || mean(B) <= 0) return(c(NA_real_, NA_real_))
  r <- mean(A) / mean(B)
  Z <- (A - r * B) / mean(B)
  c(r, if (length(A) > 1) sd(Z) / sqrt(length(A)) else NA_real_)
}

# ---- MSE ratios, bias, definedness ------------------------------------------
mse_rows <- list()
for (cu in split(units, units$cell)) {
  fid <- factor(cu$field)
  eH  <- cu$op_H - cu$h_n
  eHi <- cu$op_H - cu$h_int
  B   <- tapply(eH^2, fid, sum); Bi <- tapply(eHi^2, fid, sum)
  for (e in ests) {
    op <- cu[[paste0("op_", e)]]; raw <- cu[[paste0("est_", e)]]; d <- cu[[paste0("def_", e)]]
    eT <- op - cu$h_n
    A  <- tapply(eT^2, fid, sum)
    Ai <- tapply((op - cu$h_int)^2, fid, sum)
    Ac <- tapply(ifelse(d, (raw - cu$h_n)^2, 0), fid, sum)
    Bc <- tapply(ifelse(d, eH^2, 0), fid, sum)
    u  <- ratio_mcse(A, B); ui <- ratio_mcse(Ai, Bi)
    uc <- if (sum(Bc) > 0) ratio_mcse(Ac, Bc) else c(NA_real_, NA_real_)
    mse_rows[[length(mse_rows) + 1]] <- data.frame(
      cu[1, keys], estimator = e, F = nlevels(fid), N = nrow(cu),
      def_rate = mean(d), mean_h_n = mean(cu$h_n),
      mse = mean(eT^2), mse_H = mean(eH^2),
      ratio = u[1], mcse = u[2],
      ratio_cond = uc[1], mcse_cond = uc[2],
      ratio_int = ui[1], mcse_int = ui[2],
      bias = mean(eT), bias_cond = if (any(d)) mean(raw[d] - cu$h_n[d]) else NA_real_,
      rel_bias = mean(eT / cu$h_n),
      row.names = NULL, stringsAsFactors = FALSE)
  }
}
mse <- do.call(rbind, mse_rows)
mse <- mse[order(mse$cell, match(mse$estimator, EST_ALL)), ]
rownames(mse) <- NULL

# ---- stopping rule ------------------------------------------------------------
stop_rows <- lapply(split(mse, mse$cell), function(m) {
  m <- m[m$estimator != "H", ]
  inband <- !is.na(m$ratio) & m$ratio >= SIM$ratio_band[1] & m$ratio <= SIM$ratio_band[2]
  mx <- if (any(inband)) max(m$mcse[inband], na.rm = TRUE) else NA_real_
  Fc <- m$F[1]
  need <- any(inband) && isTRUE(mx > SIM$mcse_target) && Fc < SIM$F_max
  data.frame(m[1, keys], F = Fc, n_inband = sum(inband),
             inband_estimators = paste(m$estimator[inband], collapse = ";"),
             max_mcse_inband = mx, needs_more = need,
             next_block = if (need) Fc %/% SIM$F_block + 1L else NA_integer_,
             row.names = NULL, stringsAsFactors = FALSE)
})
stopping <- do.call(rbind, stop_rows)

# completeness: contiguous field ids 1..F, F a multiple of the block size, U units per field
comp <- do.call(rbind, lapply(split(units, units$cell), function(cu) {
  f <- sort(unique(cu$field)); Fc <- length(f)
  data.frame(cell = cu$cell[1], contiguous = identical(as.integer(f), seq_len(Fc)),
             full_blocks = Fc %% SIM$F_block == 0,
             units_ok = all(table(cu$field) == SIM$U), stringsAsFactors = FALSE)
}))
stopping <- merge(stopping, comp, by = "cell", sort = TRUE)
na_mcse <- vapply(split(mse, mse$cell), function(m) {
  m <- m[m$estimator != "H", ]
  inb <- !is.na(m$ratio) & m$ratio >= SIM$ratio_band[1] & m$ratio <= SIM$ratio_band[2]
  any(inb & is.na(m$mcse))
}, TRUE)
problems <- c(
  if (any(!stopping$contiguous)) paste("non-contiguous field ids in cells", paste(stopping$cell[!stopping$contiguous], collapse = ",")),
  if (any(!stopping$full_blocks)) paste("incomplete blocks in cells", paste(stopping$cell[!stopping$full_blocks], collapse = ",")),
  if (any(!stopping$units_ok)) paste("fields with a unit count other than U in cells", paste(stopping$cell[!stopping$units_ok], collapse = ",")),
  if (any(na_mcse)) paste("NA MCSE for an in-band estimator in cells", paste(names(na_mcse)[na_mcse], collapse = ",")))
if (length(problems)) {
  if (!partial) stop(paste(problems, collapse = "; "))
  warning(paste(problems, collapse = "; "), call. = FALSE)
}
if (!partial) {
  nm <- stopping[stopping$needs_more, ]
  writeLines(if (nrow(nm)) sprintf("%d %d", nm$cell, nm$next_block) else character(0),
             file.path(out_dir, "tasklist_next.txt"))
}

# ---- crossing n (balanced grid) ---------------------------------------------
bal <- mse[!mse$unbalanced & mse$estimator != "H", ]
cross_rows <- lapply(split(bal, list(bal$family, bal$setting, bal$tau, bal$estimator), drop = TRUE), function(m) {
  m <- m[order(m$n), ]; m <- m[is.finite(m$ratio) & m$ratio > 0, ]
  interp <- function(j) exp(log(m$n[j]) + (0 - log(m$ratio[j])) *
                              (log(m$n[j + 1]) - log(m$n[j])) / (log(m$ratio[j + 1]) - log(m$ratio[j])))
  below <- m$ratio < 1
  up <- which(below[-length(below)] & !below[-1]); dn <- which(!below[-length(below)] & below[-1])
  data.frame(family = m$family[1], setting = m$setting[1], tau = m$tau[1], estimator = m$estimator[1],
             n_values = paste(m$n, collapse = ";"),
             pattern = if (!nrow(m)) NA_character_ else if (all(below)) "below 1 at all n" else
               if (!any(below)) "at or above 1 at all n" else "crosses",
             n_cross_up = if (length(up)) interp(up[1]) else NA_real_,
             n_cross_down = if (length(dn)) interp(dn[1]) else NA_real_,
             n_sign_changes = length(up) + length(dn),
             largest_n_below_1 = if (any(below)) max(m$n[below]) else NA_real_,
             row.names = NULL, stringsAsFactors = FALSE)
})
crossing <- do.call(rbind, cross_rows)

# ---- threshold decisions -----------------------------------------------------
dec$absdist <- abs(dec$dist)
dec$bin <- cut(dec$absdist, c(0, 0.5, 1, 2, 3), include.lowest = TRUE)
dec_rate <- function(err, field) {
  s <- tapply(err, field, sum); k <- tapply(err, field, length)
  r <- sum(s) / sum(k)
  c(r, if (length(s) > 1) sd(s - r * k) / mean(k) / sqrt(length(s)) else NA_real_)
}
dests <- intersect(ests, sub("^dec_", "", grep("^dec_", names(dec), value = TRUE)))
drows <- list(); brows <- list()
for (dc in split(dec, dec$cell)) {
  for (e in dests) {
    err <- dc[[paste0("dec_", e)]] != dc$delta
    a <- dec_rate(err, dc$field)
    near <- dc$absdist <= 1
    b <- if (any(near)) dec_rate(err[near], dc$field[near]) else c(NA, NA)
    drows[[length(drows) + 1]] <- data.frame(dc[1, keys], estimator = e, n_decisions = length(err),
      err_rate = a[1], mcse = a[2], err_rate_near = b[1], mcse_near = b[2], n_near = sum(near),
      row.names = NULL, stringsAsFactors = FALSE)
    for (lv in levels(dc$bin)) {
      w <- which(dc$bin == lv)
      if (!length(w)) next
      brows[[length(brows) + 1]] <- data.frame(dc[1, keys], estimator = e, absdist_bin = lv,
        n_decisions = length(w), err_rate = mean(err[w]), row.names = NULL, stringsAsFactors = FALSE)
    }
  }
}
decision <- do.call(rbind, drows)
decision$small_n <- decision$n <= 32
decision_bin <- do.call(rbind, brows)

results <- list(mse = mse, stopping = stopping, crossing = crossing, decision = decision,
                decision_bin = decision_bin,
                meta = list(files = basename(files), n_units = nrow(units), created = Sys.time(),
                            estimators = ests, primary = EST_PRIMARY, R = R.version.string))
saveRDS(results, file.path(out_dir, "results.rds"))
write.csv(mse, file.path(out_dir, "summary.csv"), row.names = FALSE)

cat(sprintf("%d files, %d cells, %d unit rows, %d decision rows\n", length(files),
            length(unique(units$cell)), nrow(units), nrow(dec)))
cat(sprintf("stopping: %d of %d cells need another block (%s)\n", sum(stopping$needs_more), nrow(stopping),
            if (partial) "partial run: tasklist not written" else "tasklist_next.txt"))

# ---- figure -------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  pe <- setdiff(EST_PRIMARY, "H")
  pd <- mse[!mse$unbalanced & mse$estimator %in% pe & is.finite(mse$ratio), ]
  pd$estimator <- factor(pd$estimator, levels = pe)
  pd$family <- factor(pd$family, levels = SIM$families)
  pd$row <- factor(sprintf("%s, tau = %.1f", pd$setting, pd$tau),
                   levels = c(outer(SIM$settings, SIM$taus, function(s, t) sprintf("%s, tau = %.1f", s, t))))
  pal <- c(Pareto = "#2a78d6", LN = "#eb6834", Burr = "#1baf7a", NB2 = "#eda100", LNpool = "#e87ba4",
           Burrpool = "#008300", BurrPP = "#4a3aa7", HW2 = "#e34948", HW4 = "#e34948")
  lty <- c(Pareto = "solid", LN = "solid", Burr = "solid", NB2 = "solid", LNpool = "solid",
           Burrpool = "solid", BurrPP = "solid", HW2 = "solid", HW4 = "22")
  g <- ggplot(pd, aes(n, ratio, colour = estimator, linetype = estimator)) +
    geom_hline(yintercept = 1, colour = "grey45", linewidth = 0.4) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.4) +
    scale_x_log10(breaks = SIM$n_grid) + scale_y_log10() +
    scale_colour_manual(values = pal[pe], drop = FALSE) +
    scale_linetype_manual(values = lty[pe], drop = FALSE) +
    facet_grid(row ~ family, drop = FALSE) +
    labs(x = "record size n (log scale)", y = "MSE ratio vs H (operational estimator, log scale)",
         colour = NULL, linetype = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92"),
          legend.position = "top", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.background = element_rect(fill = "grey95", colour = NA))
  ggsave(fig, g, width = 14, height = 8, dpi = 150)
  cat("figure:", fig, "\n")
} else cat("ggplot2 not installed: figure skipped\n")
