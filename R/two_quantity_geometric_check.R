# ======================================================================
# two_quantity_geometric_check.R
#
# Purpose  :
#   Sufficiency check. Verifies that for geometric citations the population
#   h-index is a function of just two summaries: the number of papers n and the
#   total citation count C.
#
# Produces : console check (no figure)
# Reads    : nothing (self-contained simulation)
# Requires : base R only
# Run      : from the supplement root --  Rscript two_quantity_geometric_check.R
# Seed     : set.seed(1729) (all stochastic scripts use this seed)
# ======================================================================

# Geometric citations X ~ Geom(q) on {0,1,2,...}: P(X=x)=q(1-q)^x, S(h)=(1-q)^h.
# Claim: population h_n solves n (1-q)^h = h, i.e. h_n = W(a n)/a, a = -log(1-q).
# And a is a function of just (n, C): q_hat = n/(n+C), 1-q_hat = C/(n+C),
#   a_hat = log(1 + n/C).  => h is a function of TWO productivity quantities (n, C).
lambertW0 <- function(z){ w <- log(z+1); for(i in 1:60){ e<-exp(w); w <- w - (w*e - z)/(e*(w+1)) }; w }
set.seed(1729)
chk <- function(n,q){
  a <- -log(1-q); hpop <- lambertW0(a*n)/a
  X <- matrix(rgeom(2000*n,q),nrow=2000)
  H <- apply(X,1,function(x){xs<-sort(x,decreasing=TRUE);m<-which(xs>=seq_along(xs));if(length(m))max(m) else 0})
  C <- mean(rowSums(X)); a_hat <- log(1+n/C); hpop_2q <- lambertW0(a_hat*n)/a_hat
  cat(sprintf("n=%4d q=%.2f | h_pop(q)=%.2f  E[H_n]=%.2f  h_2q(n,Cbar)=%.2f\n",
              n,q,hpop,mean(H),hpop_2q))
}
for(q in c(0.1,0.3)) for(n in c(50,200,1000)) chk(n,q)
