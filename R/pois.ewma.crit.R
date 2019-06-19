# Computation of Poisson EWMA control limits
pois.ewma.crit <- function(lambda, L0, mu0, z0, AU=3, sided="two", design="sym", rando=FALSE, mcdesign="transfer", N=101, jmax=4) {
  if ( lambda <= 0 | lambda > 1 )       stop("lambda has to be between 0 and 1")
  if ( L0 < 1 )                         stop("L0 has to be larger")
  if ( mu0 < 0 )     	                stop("wrong value for mu0")
  ctyp <- pmatch(sided, c("upper", "lower", "two")) - 1
  if ( is.na(ctyp) )                    stop("invalid ewma type")
  dtyp <- pmatch(design, c("sym", "unb")) - 1
  if ( is.na(dtyp) )                    stop("invalid design type")
  mcd <- pmatch(mcdesign, c("classic", "transfer")) - 1
  if ( is.na(mcd) )                    stop("invalid mcdesign value")
  LL <- ifelse(dtyp==0, 1, 2)
  if (dtyp==1 & rando) LL <- 4
  crit <- .C("cewma_crit_be",
             as.integer(ctyp), as.integer(dtyp), as.integer(mcd), as.integer(rando), as.double(lambda), as.double(L0), as.double(AU),
             as.double(mu0), as.double(z0), as.integer(N), as.integer(jmax),
             ans=double(length=LL), PACKAGE="spc")$ans
  if (dtyp==1 & rando) {
    names(crit) <- c("AL", "AU", "gL", "gU")
  } else {
    if ( ctyp==0 ) names(crit) <- "AU"
    if ( ctyp==1 ) names(crit) <- "AL"
    if ( dtyp==0 ) names(crit) <- "A"
    if ( dtyp==1 ) names(crit) <- c("AL", "AU")
  }
  crit
}
