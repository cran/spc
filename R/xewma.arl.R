# Computation of EWMA ARLs (mean monitoring)
xewma.arl <- function(l, cE, mu, zr=0, hs=0, sided="one", limits="fix", q=1, steady.state.mode="conditional", r=40) {
  if ( l<=0 | l>2 )
    stop("l has to be between 0 and 2")
  if ( any(cE<=0) )
    warning("usually, cE has to be positive")
  if ( limits!="cfar" ) {  
    if ( zr>cE & sided=="one" )
      stop("wrong reflexion border")
    if ( (sided=="two" & abs(hs)>cE) | (sided=="one" & (hs<zr | hs>cE)) ) 
      warning("unusual headstart")
  }
  if ( r<4 )
    stop("r is too small")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if ( is.na(ctyp) )
    stop("invalid ewma type")
  ltyp <- -1 + pmatch(limits,
          c("fix", "vacl", "fir", "both", "Steiner", "stat", "cfar", "limit", "fixW", "fixC"))            
  if ( is.na(ltyp) )
    stop("invalid limits type")
  if ( (sided=="one") & !(limits %in% c("fix", "vacl", "stat", "limit", "fixW")) )
    stop("not supported for one-sided EWMA (not reasonable or not implemented yet")
  q <- round(q)
  if ( q<1 )
    stop("wrong change point position (q)")
  styp <- pmatch(steady.state.mode, c("conditional", "cyclical")) - 1
  if (is.na(styp))
    stop("invalid steady.state.mode")
  if ( limits=="cfar" ) {
    ctyp <- length(cE) # nc
    hs <- cE[length(cE)] # cinf
  }
  if ( limits=="fix" & q>1 & styp==0 ) {
    arl <- .C("xewma_arl",as.integer(ctyp),as.double(l),
              as.double(cE),as.double(zr),as.double(hs),
              as.double(mu),as.integer(ltyp),as.integer(r),as.integer(q),as.integer(styp),
              ans=double(length=q), PACKAGE="spc")$ans 
  } else {
    arl <- .C("xewma_arl",as.integer(ctyp),as.double(l),
              as.double(cE),as.double(zr),as.double(hs),
              as.double(mu),as.integer(ltyp),as.integer(r),as.integer(q),as.integer(styp),
              ans=double(length=1), PACKAGE="spc")$ans
  }
  names(arl) <- NULL
  return (arl)
}
