# Computation of Poisson CUSUM ARLs
pois.cusum.arl <- function(mu, km, hm, m, i0=0, sided="upper", rando=FALSE, gamma=0, km2=0, hm2=0, m2=0, i02=0, gamma2=0) {
  if ( mu < 0 )         stop("mu has to be positive")
  if ( km < 1 )         stop("km has to be >= 1") 
  if ( hm < 1 )         stop("hm has to be >= 1")
  if (  m < 1 )         stop("m has to be >= 1")
  if ( !(i0 %in% 0:hm) )     stop("head start i0 has to be an integer in 0, 1, ..., hm")
  ctyp <- pmatch(sided, c("upper", "lower", "two")) - 1
  if ( is.na(ctyp) )    stop("invalid cusum type")
  if ( rando ) {
    if ( gamma < 0 | gamma > 1)   stop("gamma has to be a probability value (0<=.<=1)")
  }
  if ( ctyp==2 ) {
    if ( km2 < 1 )         stop("km2 has to be >= 1") 
    if ( hm2 < 1 )         stop("hm2 has to be >= 1")
    if (  m2 < 1 )         stop("m2 has to be >= 1")
    if ( !(i02 %in% 0:hm) )     stop("head start i02 has to be an integer in 0, 1, ..., hm")
    if ( rando ) {
      if ( gamma2 < 0 | gamma2 > 1)   stop("gamma2 has to be a probability value (0<=.<=1)")
    }
  }
  arl <- .C("ccusum_arl_be",
            as.integer(ctyp), as.integer(rando), as.double(mu),
            as.integer(km), as.integer(hm), as.integer(m), as.integer(i0), as.double(gamma),
            as.integer(km2), as.integer(hm2), as.integer(m2), as.integer(i02), as.double(gamma2),
            ans=double(length=1), PACKAGE="spc")$ans 
  names(arl) <- "arl"
  arl
}
