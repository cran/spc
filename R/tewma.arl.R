# Computation of TEWMA ARLs for iid Poisson
tewma.arl <- function(lambda, k, lk, uk, mu, z0, rando=FALSE, gl=0, gu=0) {
  if ( lambda <= 0 || lambda > 1 )      stop("lambda has to be between 0 and 1")
  if ( k < 1 )		                    stop("k must be >= 1")
  if ( mu < 0 )		                    stop("mu must be positive")  
  if ( lk > z0 | z0 > uk )              stop("wrong headstart")  
  if ( 0 > gl | gl > 1 )                stop("wrong value for gl")
  if ( 0 > gu | gu > 1 )                stop("wrong value for gu")
  irando <- as.numeric(rando)
  
  arl <- .C("tewma_arl_wowR",
            as.integer(irando), as.double(lambda),
            as.integer(k), as.integer(lk), as.integer(uk), as.double(gl), as.double(gu),
            as.double(z0), as.double(mu), ans=double(length=1), PACKAGE="spc")$ans 
  names(arl) <- "arl"
  arl
}
