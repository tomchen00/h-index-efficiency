# ======================================================================
# sim_decisions.R
#
# Purpose  : Threshold decisions at short record sizes, conditional on the observed margin H - t
#            (as in the ASN application, where applicants are selected on their observed index).
#            For balanced cells with n <= 32 and each estimator T: P(delta = 1 | margin), the error
#            rates of the rules 1{H >= t} and 1{T >= t}, and, among rows where the two rules disagree,
#            the share in which T agrees with the true decision delta = 1{h_n >= t}.
# Reads    : data/sim/raw/*.rds, data/sim/cells.csv
# Produces : data/sim/decisions_short_n.rds, data/sim/decisions_short_n.csv
# Run      : from the repository root -- Rscript R/sim_decisions.R
# ======================================================================
ests <- c("LN", "NB2", "LNpool", "Burrpool", "BurrPP")
files <- list.files(file.path("data", "sim", "raw"), pattern = "\\.rds$", full.names = TRUE)
cells <- read.csv(file.path("data", "sim", "cells.csv"))
keep_cells <- cells$cell[cells$n <= 32 & !cells$unbalanced]
files <- files[as.integer(sub("cell_([0-9]+)_block.*", "\\1", basename(files))) %in% keep_cells]
acc <- list()
for (f in files) {
  x <- readRDS(f)
  u <- x$units[, c("field", "unit", "est_H")]
  d <- merge(x$dec, u, by = c("field", "unit"))
  d$margin <- pmax(-2, pmin(2, d$est_H - d$t))
  for (e in ests) {
    de <- d[[paste0("dec_", e)]]; dh <- d$dec_H
    key <- interaction(d$margin, drop = TRUE)
    g <- split(data.frame(delta = d$delta, dh = dh, de = de), d$margin)
    for (m in names(g)) {
      G <- g[[m]]
      mc <- x$meta$cell; id <- paste(mc$family, mc$setting, mc$tau, mc$n, e, m, sep = "|")
      a <- acc[[id]]; if (is.null(a)) a <- c(rows = 0, pass = 0, errH = 0, errT = 0, disagree = 0, T_right = 0)
      dis <- G$dh != G$de
      a <- a + c(nrow(G), sum(G$delta), sum(G$dh != G$delta), sum(G$de != G$delta), sum(dis), sum(dis & G$de == G$delta))
      acc[[id]] <- a
    }
  }
}
R <- do.call(rbind, lapply(names(acc), function(k) {
  p <- strsplit(k, "\\|")[[1]]; a <- acc[[k]]
  data.frame(family = p[1], setting = p[2], tau = as.numeric(p[3]), n = as.integer(p[4]), estimator = p[5],
             margin = as.integer(p[6]), rows = a[["rows"]], p_pass = a[["pass"]] / a[["rows"]],
             err_H = a[["errH"]] / a[["rows"]], err_T = a[["errT"]] / a[["rows"]],
             disagree = a[["disagree"]] / a[["rows"]],
             T_right_when_disagree = if (a[["disagree"]] > 0) a[["T_right"]] / a[["disagree"]] else NA_real_)
}))
R <- R[order(R$setting, R$family, R$tau, R$n, R$estimator, R$margin), ]
saveRDS(R, file.path("data", "sim", "decisions_short_n.rds"))
write.csv(R, file.path("data", "sim", "decisions_short_n.csv"), row.names = FALSE)
cat("rows", nrow(R), "\n")
