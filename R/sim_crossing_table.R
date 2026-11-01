# ======================================================================
# sim_crossing_table.R
#
# Purpose  :
#   Table 2 of the paper: record size n at which the mean squared error
#   ratio MSE(T)/MSE(H_n) of an estimator first rises above one, by
#   citation distribution, in the short setting on the balanced grid.
#   Uses the operational ratio (undefined fits replaced by H_n) and the
#   log-log interpolation of sim_aggregate.R between design sizes.
#   Codes: "--" = ratio below one at every simulated size,
#          "*"  = ratio at or above one at every simulated size.
#
# Produces : data/sim/crossing_table.txt (LaTeX rows); prints them
# Reads    : data/sim/summary.csv
# Requires : base R
# Run      : from the supplement root --  Rscript R/sim_crossing_table.R
# ======================================================================

s <- read.csv("data/sim/summary.csv", stringsAsFactors = FALSE)
s <- s[!s$unbalanced & s$setting == "short" & s$estimator != "H", ]
fam <- c("C1", "C2", "T1", "T2", "T3", "T4", "T5", "T7")

cross_n <- function(n, r) {
  o <- order(n); n <- n[o]; r <- r[o]; ok <- is.finite(r) & r > 0; n <- n[ok]; r <- r[ok]
  below <- r < 1
  if (all(below)) return("--")
  if (!any(below)) return("$*$")
  up <- which(below[-length(below)] & !below[-1])
  if (!length(up)) return("--")
  j <- up[1]
  x <- exp(log(n[j]) + (0 - log(r[j])) * (log(n[j + 1]) - log(n[j])) / (log(r[j + 1]) - log(r[j])))
  sprintf("%d", as.integer(floor(x + 0.5)))
}

rows <- list(c("pooled lognormal", "LNpool", 0.3), c("", "LNpool", 0.1),
             c("pooled Burr XII", "Burrpool", 0.3), c("", "Burrpool", 0.1),
             c("lognormal", "LN", 0.3), c("NB2", "NB2", 0.3))
out <- character(0)
for (rw in rows) {
  d <- s[s$estimator == rw[2] & s$tau == as.numeric(rw[3]), ]
  v <- sapply(fam, function(f) { m <- d[d$family == f, ]; cross_n(m$n, m$ratio) })
  out <- c(out, sprintf("%s & %s & %s\\\\", rw[1], rw[3], paste(v, collapse = " & ")))
}
writeLines(out, "data/sim/crossing_table.txt")
writeLines(out)
