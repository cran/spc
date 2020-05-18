# Computation of Poisson EWMA steady-state ARLs
pois.ewma.ad <- function(lambda, AL, AU, mu0, mu, sided="two", rando=FALSE, gL=0, gU=0, mcdesign="classic", N=101) {
  if ( lambda <= 0 | lambda > 1 )       stop("lambda has to be between 0 and 1")
  if ( AL < 0 | AU < 0 )                stop("control limit factors must be positive")
  if ( mu0 < 0 )     	                stop("wrong value for mu0")
  if ( mu < 0 )     	                stop("wrong value for mu")
  hL  <- mu0 - AL*sqrt(lambda*mu0/(2-lambda))
  hU  <- mu0 + AU*sqrt(lambda*mu0/(2-lambda))
  ctyp <- pmatch(sided, c("upper", "lower", "two")) - 1
  if ( is.na(ctyp) )                    stop("invalid ewma type")
  mcd <- pmatch(mcdesign, c("classic", "transfer")) - 1
  if ( is.na(mcd) )                    stop("invalid mcdesign value")
  if ( rando ) {
    if ( gL < 0 | gL > 1 )              stop("wrong value for gL") 
    if ( gU < 0 | gU > 1 )              stop("wrong value for gU")
  }
  ad <- .C("cewma_ad_be",
           as.integer(ctyp), as.integer(mcd), as.integer(rando), as.double(lambda),
           as.double(AL), as.double(AU), as.double(gL), as.double(gU),
           as.double(mu0), as.double(mu), as.integer(N),
           ans=double(length=1), PACKAGE="spc")$ans 
  names(ad) <- "ad"
  ad
}
