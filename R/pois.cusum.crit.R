# Computation of Poisson CUSUM alarm threshold and randomization constant
pois.cusum.crit <- function(mu0, km, A, m, i0=0, sided="upper", rando=FALSE) {
  if ( mu0 < 0 )        stop("mu0 has to be positive")
  if ( km < 1 )         stop("km has to be >= 1") 
  if ( A < 1 )          stop("A has to be >= 1")
  if ( m < 1 )          stop("m has to be >= 1")
  i0 <- round(i0)
  if ( i0 < 0 | i0 > 100 )     stop("head start i0 has to be an integer in 0, 1, ..., 100")
  ctyp <- pmatch(sided, c("upper", "lower", "two")) - 1
  if ( is.na(ctyp) )    stop("invalid cusum type")
  cv <- .C("ccusum_crit_be",
           as.integer(ctyp), as.integer(rando),
           as.double(mu0), as.integer(km), as.double(A), as.integer(m), as.integer(i0),
           ans=double(length=2), PACKAGE="spc")$ans 
  names(cv) <- c("hm", "gamma")
  cv[1] <- as.integer(cv[1])
  cv
}
