# Computation of I-MR ARLs
imr.arl <- function(M, Ru, mu, sigma, vsided="upper", Rl=0, cmode="coll", N=30, qm=30) {
  if ( M <= 0 )         stop("M has to be positive")
  if ( Ru <= 0 )        stop("Ru has to be positive") 
  if ( sigma <= 0 )     stop("sigma has to be positive")
  if ( Rl < 0 )         stop("Rl has to be non-negative") 
  if ( N < 2 )          stop("N has to be >= 2")
  if ( qm < 10 )        stop("qm has to be >= 10")
  vsided <- tolower(vsided)
  cmode <- tolower(cmode)  
  vtyp <- -1 + pmatch(vsided, c("upper", "two"))
  
  if ( Ru >= 2*M & vsided == "upper" ) {
    Lu <- 1 / ( 2*pnorm(-M, mean=mu, sd=sigma) )
    arl <- 1 + Lu * ( 2*pnorm(M, mean=mu, sd=sigma) - 1 )
  } else {
    arl <- -1
    if ( cmode == "coll" | vsided=="two" ) {
      arl <- .C("imr_arl", as.double(M), as.double(Rl), as.double(Ru), as.double(mu), as.double(sigma),
                           as.integer(vtyp), as.integer(N), as.integer(qm), ans=double(length=1), PACKAGE="spc")$ans 
    }
    if ( cmode == "crowder" ) {
      if ( vsided == "two" ) warning("confirmed only for upper MR")
      arl <- imr.arl.C1987b(M, Ru, mu/sigma, N=N)
    }
    if ( cmode %in% c("gl", "rectangular", "trapezoid", "simpson", "simpson3_8") ) {
      if ( vsided == "two" ) warning("confirmed only for upper MR")
      arl <- imr.arl.Ny(M, Ru, mu/sigma, N, kind=cmode)
    }
  }  
  names(arl) <- "arl"
  arl
}
