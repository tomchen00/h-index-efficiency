# ======================================================================
# sim_run.R
#
# Purpose  :
#   Runs one cell of the simulation grid for one block of fields. Each field
#   holds U = 30 units drawn around the cell's calibrated centre; for every
#   unit it records the law, the targets (h_n from the log-linear extension
#   of S^-, h_int), the empirical H and every estimator (Pareto, lognormal,
#   Burr, NB2, pooled lognormal, pooled Burr, partially pooled Burr,
#   Hill-Weissman kappa = 2, 4, untruncated lognormal/Burr sensitivity):
#   raw value (NA if undefined), definedness, operational value (fallback
#   to H), fit diagnostics, and threshold-decision rows for
#   t in {h_int - 2, ..., h_int + 2}, t >= 1.
#
# Produces : data/sim/raw/cell_<i>_block_<b>.rds  (list: units, dec, seeds, meta);
#            meta holds the fingerprint (md5 of sim_common.R, sim_run.R,
#            data/sim/cells.csv), the cell parameters, field ids and
#            estimator set. If the output exists and matches all four, the
#            run exits 0 without recomputing; if it exists and differs, it
#            stops with an error. Written to a PID-specific temporary file,
#            then renamed.
# Reads    : R/sim_common.R, data/sim/cells.csv
# Requires : base R (stats, parallel)
# Run      : from the repository root --
#              Rscript R/sim_run.R <cell_index> <field_block> [n_fields] [out_dir]
#            block b covers fields (b-1)*150 + 1, ..., (b-1)*150 + n_fields
#            (n_fields defaults to 150; out_dir defaults to data/sim/raw)
# Seed     : field f of cell i uses L'Ecuyer-CMRG stream J = (i-1)*600 + f,
#            reached from set.seed(20260930) by J calls of nextRNGStream;
#            each field's .Random.seed is stored in the output
# ======================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript R/sim_run.R <cell_index> <field_block> [n_fields] [out_dir]")
cell_i  <- as.integer(args[1])
block   <- as.integer(args[2])
source("R/sim_common.R")
n_fields <- if (length(args) >= 3) as.integer(args[3]) else SIM$F_block
out_dir  <- if (length(args) >= 4) args[4] else "data/sim/raw"
stopifnot(block >= 1, block * SIM$F_block <= SIM$F_max, n_fields >= 1, n_fields <= SIM$F_block)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cells <- read.csv("data/sim/cells.csv", stringsAsFactors = FALSE)
cell <- cells[cells$cell == cell_i, ]
if (nrow(cell) != 1) stop("cell index not found: ", cell_i)

fields <- (block - 1) * SIM$F_block + seq_len(n_fields)
fp <- sim_fingerprint()
fn <- file.path(out_dir, sprintf("cell_%d_block_%d.rds", cell_i, block))

same_cell <- function(a, b) {
  nm <- names(b)
  identical(sort(names(a)), sort(nm)) &&
    all(vapply(nm, function(v) isTRUE(all.equal(a[[v]], b[[v]], check.attributes = FALSE)), TRUE))
}
if (file.exists(fn)) {
  old <- tryCatch(readRDS(fn), error = function(e) NULL)
  if (is.null(old)) stop("existing output unreadable: ", fn)
  m <- old$meta
  ok <- identical(m$fingerprint, fp) && identical(as.integer(m$fields), as.integer(fields)) &&
    identical(m$estimators, EST_ALL) && same_cell(as.list(m$cell), as.list(cell))
  if (!ok) stop("existing output does not match the current code/cell/fields/estimators: ", fn,
                "\n  stored fingerprint: ", format(m$fingerprint), "\n  current fingerprint: ", fp)
  cat(sprintf("cell %d block %d: output exists and matches; skipped (%s)\n", cell_i, block, fn))
  quit(save = "no", status = 0)
}

seed <- sim_stream_seed(cell_i, fields[1])        # stream of the first field in the block

units <- vector("list", n_fields); dec <- vector("list", n_fields)
seeds <- matrix(NA_integer_, n_fields, length(seed))
secs  <- numeric(n_fields)
t_start <- Sys.time()
for (j in seq_len(n_fields)) {
  if (j > 1) seed <- parallel::nextRNGStream(seed)
  seeds[j, ] <- seed
  assign(".Random.seed", seed, envir = .GlobalEnv)
  t0 <- proc.time()[3]
  res <- sim_field(cell)
  secs[j] <- proc.time()[3] - t0
  key <- data.frame(cell = cell_i, family = cell$family, setting = cell$setting, tau = cell$tau,
                    n = cell$n, unbalanced = as.logical(cell$unbalanced), block = block,
                    field = fields[j], stringsAsFactors = FALSE)
  units[[j]] <- cbind(key[rep(1, nrow(res$units)), ], res$units, row.names = NULL)
  if (!is.null(res$dec))
    dec[[j]] <- cbind(key[rep(1, nrow(res$dec)), ], res$dec, row.names = NULL)
}

out <- list(
  units = do.call(rbind, units),
  dec   = do.call(rbind, dec),
  seeds = data.frame(field = fields, seeds),
  meta  = list(fingerprint = fp, cell = cell, block = block, fields = fields, U = SIM$U,
               estimators = EST_ALL, primary = EST_PRIMARY,
               sec_per_field = secs, wall_start = t_start, wall_end = Sys.time(),
               R = R.version.string, rng = RNGkind(), master_seed = SIM$master_seed)
)
tmp <- sprintf("%s.tmp%d", fn, Sys.getpid())
saveRDS(out, tmp)
if (!file.rename(tmp, fn)) stop("could not write ", fn)
cat(sprintf("cell %d (%s %s tau=%.1f n=%d%s) block %d: %d fields, %.1f s/field, total %.1f s -> %s\n",
            cell_i, cell$family, cell$setting, cell$tau, cell$n,
            if (isTRUE(as.logical(cell$unbalanced))) " unbal" else "",
            block, n_fields, mean(secs), sum(secs), fn))
