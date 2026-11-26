# asn_calibration_targets.R: the two calibration targets of the simulation study (Section 5.1, Supplementary Table S5).
# Median over the habilitation records (prorata version) of H_n/n and of the ratio q90/q50 of sample quantiles,
# with R's default quantile (type 7, linear interpolation between order statistics); records whose sample median is 0
# give no finite ratio and are excluded. Output (2026-11-25): median H_n/n = 0.444, median q90/q50 = 3.756 (115 excluded).
# Run from the repository root: Rscript R/asn_calibration_targets.R
suppressPackageStartupMessages(library(jsonlite))
A <- readRDS("data/asn/derived/asn_estimates_prorata.rds")
recs <- lapply(readLines("data/asn/derived/asn_records.jsonl"), fromJSON)
ids <- vapply(recs, `[[`, "", "id")
X <- lapply(recs, function(r) floor(as.numeric(r$x_prorata)))[match(A$id, ids)]
r <- vapply(X, function(x) { q <- quantile(x, c(0.5, 0.9), type = 7); q[[2]] / q[[1]] }, 0)
r[!is.finite(r)] <- NA
cat(sprintf("records %d; median n %g; median H_n/n %.3f; median q90/q50 %.3f (%d excluded)\n",
            nrow(A), median(A$n), median(A$H / A$n), median(r, na.rm = TRUE), sum(is.na(r))))
