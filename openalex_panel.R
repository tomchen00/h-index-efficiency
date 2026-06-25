# ======================================================================
# openalex_panel.R
#
# Purpose  :
#   OpenAlex data pipeline, step 1. Fetches the six long-career authors' full
#   citation records, then computes each author's empirical H, Pareto plug-in,
#   winsorized plug-in, and bootstrap variance ratios. The raw pull is cached so
#   the rest runs offline.
#
# Produces : data/authors_data.rds (raw pull), data/panel_summary.rds; console table
# Reads    : data/authors_data.rds if present (uses the cache instead of re-fetching)
# Requires : jsonlite; internet only on the first run (the shipped cache avoids it)  (CRAN packages)
# Run      : from the supplement root --  Rscript openalex_panel.R
# Seed     : set.seed(1729) (all stochastic scripts use this seed)
# ======================================================================

## ensure output folders exist (some zip extractors drop empty dirs)
for (d_ in c("data","figures")) if (!dir.exists(d_)) dir.create(d_, showWarnings=FALSE)

# =====================================================================
# Multi-author OpenAlex panel for the JOI illustration.
# For a cross-field set of long-career authors:
#   - pull works (cited_by_count, year); CACHE to authors_data.rds
#   - empirical H, Pareto-MLE plug-in h_hat, winsorized h_c
#   - tail goodness-of-fit: Hill-plot stability over x_min + KS distance
#   - bootstrap Var(H), Var(h_hat), Var(h_c) -> efficiency gap per author
# Self-contained (no poweRlaw). Run once; the .rds caches the pulled data.
# =====================================================================
suppressMessages({library(jsonlite)})
set.seed(1729)
DATADIR <- "data"
RDS <- file.path(DATADIR, "authors_data.rds")

resolve_author <- function(name){
  u <- sprintf("https://api.openalex.org/authors?search=%s&per-page=10&mailto=anonymous@example.org",
               URLencode(name, reserved=TRUE))
  j <- tryCatch(fromJSON(u), error=function(e) NULL)
  if (is.null(j) || length(j$results)==0) return(NULL)
  r <- j$results
  i <- which.max(r$works_count)              # most prolific match = the famous one
  list(id=sub(".*/", "", r$id[i]), name=r$display_name[i],
       works=r$works_count[i], cites=r$cited_by_count[i])
}
fetch_works <- function(author_id, max_n=6000){
  base <- sprintf("https://api.openalex.org/works?filter=author.id:%s&per-page=200&select=cited_by_count,publication_year&cursor=", author_id)
  cur <- "*"; cc<-c(); yr<-c()
  repeat{
    u <- paste0(base, URLencode(cur, reserved=TRUE), "&mailto=anonymous@example.org")
    j <- tryCatch(fromJSON(u), error=function(e) NULL)
    if (is.null(j) || length(j$results)==0) break
    cc<-c(cc,j$results$cited_by_count); yr<-c(yr,j$results$publication_year)
    cur<-j$meta$next_cursor; if(is.null(cur)||length(cc)>=max_n) break
  }
  data.frame(cites=as.integer(cc), year=as.integer(yr))
}

panel <- list(
  Economics   = "Daron Acemoglu",
  Physics     = "Edward Witten",
  Mathematics = "Terence Tao",
  CS_ML       = "Geoffrey Hinton",
  Statistics  = "Robert Tibshirani",
  Genomics    = "Eric Lander"
)

if (file.exists(RDS)) {
  cat("Loading cached", RDS, "\n"); A <- readRDS(RDS)
} else {
  A <- list()
  for (fld in names(panel)){
    info <- resolve_author(panel[[fld]])
    if (is.null(info)){ cat("FAILED resolve:", panel[[fld]],"\n"); next }
    Sys.sleep(0.3)
    w <- fetch_works(info$id)
    cat(sprintf("%-12s %-22s id=%s works=%d pulled=%d maxcite=%d\n",
                fld, info$name, info$id, info$works, nrow(w), max(w$cites,na.rm=TRUE)))
    A[[fld]] <- list(field=fld, info=info, works=w)
    Sys.sleep(0.3)
  }
  saveRDS(A, RDS); cat("Cached ->", RDS, "\n")
}

# ---- estimators + diagnostics ----
hindex <- function(x){ xs<-sort(x,decreasing=TRUE); w<-which(xs>=seq_along(xs)); if(length(w)) max(w) else 0 }
pareto_a <- function(xt, xmin=1){ z<-xt[xt>=xmin]; length(z)/sum(log(z/xmin)) }   # known-scale-ish MLE
hplug <- function(a,m) m^(1/(1+a))
g_win <- function(a,c) (1/a)*(1-exp(-a*c))
alpha_win <- function(logx,c){ Wb<-mean(pmin(logx,c)); if(Wb<=0||Wb>=c) return(length(logx)/sum(logx))
  uniroot(function(a) g_win(a,c)-Wb, c(1e-3,50))$root }
# KS distance between empirical tail (x>=xmin) and fitted Pareto
ks_pareto <- function(xt, xmin=1){ z<-sort(xt[xt>=xmin]); n<-length(z); a<-pareto_a(z,xmin)
  Fhat<-1-(z/xmin)^(-a); Femp<-(seq_len(n))/n; max(abs(Femp-Fhat)) }

cat("\n", strrep("=",92), "\n", sep="")
cat(sprintf("%-12s %5s %5s %6s %6s %7s %7s | boot Var ratios: %7s %7s\n",
            "field","n","H","alpha","KS",".hhat",".hc","hhat/H","hc/H"))
summ <- data.frame()
for (fld in names(A)){
  w <- A[[fld]]$works; x <- w$cites[!is.na(w$cites)]
  xt <- x[x>=1]; m <- length(xt)
  H <- hindex(x); a <- pareto_a(xt); hh <- hplug(a,m)
  cq <- quantile(log(xt),0.95); hc <- hplug(alpha_win(log(xt),cq), m)
  ks <- ks_pareto(xt)
  # bootstrap
  B<-8000; bH<-bHh<-bHc<-numeric(B); zeros<-x[x<1]
  for (b in 1:B){ xb<-sample(xt,m,replace=TRUE)
    bH[b]<-hindex(c(xb,zeros)); ab<-pareto_a(xb); bHh[b]<-hplug(ab,m)
    cqb<-quantile(log(xb),0.95); bHc[b]<-hplug(alpha_win(log(xb),cqb),m) }
  rH<-var(bH); rHh<-var(bHh); rHc<-var(bHc)
  cat(sprintf("%-12s %5d %5d %6.3f %6.3f %7.1f %7.1f | %7.3f %7.3f\n",
              fld, m, H, a, ks, hh, hc, rHh/rH, rHc/rH))
  summ <- rbind(summ, data.frame(field=fld, author=A[[fld]]$info$name, n=m, H=H,
                alpha=a, KS=ks, hhat=hh, hc=hc, varH=rH, varhh=rHh, varhc=rHc,
                ratio_hh=rHh/rH, ratio_hc=rHc/rH))
}
cat(strrep("=",92),"\n")
saveRDS(summ, file.path(DATADIR,"panel_summary.rds"))
cat("\nHill-plot stability (alpha_hat vs x_min) for first author:\n")
w1<-A[[1]]$works; x1<-w1$cites[w1$cites>=1]
for (xm in c(1,2,5,10,20,50)) cat(sprintf("  x_min=%3d  n_tail=%5d  alpha=%.3f\n", xm, sum(x1>=xm), pareto_a(x1,xm)))
cat("\nDONE.\n")
