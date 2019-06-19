# Computation of 'Euklid'-EWMA ARLs for iid Poisson (Rakitzis/Castagliola/Maravelakis 2015)
euklid.ewma.arl <- function(gX, gY, kL, kU, mu, y0, r0=0) {
  if ( gX <= 0 )            stop("gX must be positive")
  if ( gY <= 0 )            stop("gY must be positive")
  if ( kL <= 0 )            stop("kL must be positive")
  if ( kU <= 0 )            stop("kU must be positive")
  if ( kU < kL )            stop("kU must be larger than or equal to kL")
  if ( mu <= 0 )            stop("mu must be positive")  
  if ( y0 <= 0 )            stop("y0 must be positive")
  if ( r0 < 0 )             stop("r0 must be non-negative")

  arl <- .C("euklid_ewma_arl",
            as.integer(gX), as.integer(gY), as.integer(kL), as.integer(kU),
            as.double(mu), as.double(y0), as.integer(r0),
            ans=double(length=1), PACKAGE="spc")$ans 
  names(arl) <- "arl"
  arl
}
