# Computation of attribute p EWMA ARLs
p.ewma.arl <- function(lambda, ucl, n, p, z0, sided="upper", lcl=NULL, d.res=1, r.mode="ieee.round", i.mode="integer") {
  i.r.mode <- -2 + pmatch(r.mode, c("gan.floor", "floor", "ceil", "ieee.round", "round", "mix"))
  i.i.mode <- -1 + pmatch(i.mode, c("integer", "half"))
  if ( lambda <= 0 || lambda > 1 )      stop("lambda has to be between 0 and 1")
  if ( ucl < 0 )                        stop("ucl must be larger than 0")
  if ( sided == "upper" ) lcl <- 0
  if ( sided == "lower" ) ucl <- n
  if ( sided == "two" ) {
    if ( is.null(lcl) )                 stop("lcl must be set")
    if ( lcl < 0 )                      stop("lcl must be larger than 0")
  }
  if ( n < 1 )		                    stop("n must be >= 1")
  if ( 0 > p | p > 1 )	                stop("wrong value for p")
  if ( is.null(lcl) ) lcl <- 0
  if ( z0 < lcl | z0 > ucl )            stop("wrong headstart")
  ctyp <- pmatch(sided, c("upper", "lower", "two")) - 1
  if ( is.na(ctyp) )                    stop("invalid ewma type")
  if ( d.res < 1 )                      stop("d.res too small")
  if ( is.na(i.r.mode) )                stop("invalid round mode")
  if ( is.na(i.i.mode) )                stop("invalid interval mode")
  arl <- .C("ewma_p_arl_be",
            as.integer(ctyp), as.double(lambda), as.double(lcl), as.double(ucl), as.integer(n),
            as.double(p), as.double(z0), as.integer(d.res), as.integer(i.r.mode), as.integer(i.i.mode),
            ans=double(length=1), PACKAGE="spc")$ans 
  names(arl) <- "arl"
  arl
}
