# Computation of Poisson EWMA ARLs
pois.ewma.arl <- function(lambda, AL, AU, mu0, z0, mu, sided="two", rando=FALSE, gL=0, gU=0, mcdesign="transfer", N=101) {
  if ( lambda <= 0 | lambda > 1 )       stop("lambda has to be between 0 and 1")
  if ( AL < 0 | AU < 0 )                stop("control limit factors must be positive")
  if ( mu0 < 0 )     	                stop("wrong value for mu0")
  if ( mu < 0 )     	                stop("wrong value for mu")
  hL  <- mu0 - AL*sqrt(lambda*mu0/(2-lambda))
  hU  <- mu0 + AU*sqrt(lambda*mu0/(2-lambda))
  if ( z0 < hL | z0 > hU )             stop("wrong headstart")
  ctyp <- pmatch(sided, c("upper", "lower", "two")) - 1
  if ( is.na(ctyp) )                    stop("invalid ewma type")
  mcd <- pmatch(mcdesign, c("classic", "transfer")) - 1
  if ( is.na(mcd) )                    stop("invalid mcdesign value")
  if ( rando ) {
    if ( gL < 0 | gL > 1 )              stop("wrong value for gL") 
    if ( gU < 0 | gU > 1 )              stop("wrong value for gU")
  }
  arl <- .C("cewma_arl_be",
            as.integer(ctyp), as.integer(mcd), as.integer(rando), as.double(lambda),
            as.double(AL), as.double(AU), as.double(gL), as.double(gU),
            as.double(mu0), as.double(z0), as.double(mu), as.integer(N),
            ans=double(length=1), PACKAGE="spc")$ans 
  names(arl) <- "arl"
  arl
}
