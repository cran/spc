# Computation of EWMA critical values for given ARL (mean monitoring)
xewma.crit <- function(l,L0,mu0=0,zr=0,hs=0,sided="one",limits="fix",r=40,c0=NULL,nmax=1e4) {
  if ( l<=0 | l>2 ) 
    stop("l has to be between 0 and 2")
  if ( L0<1 ) 
    stop("L0 is too small")
  if ( r<4 )
    stop("r is too small")
  if ( sided=="one" & hs<zr )
    warning("unusual headstart")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if ( is.na(ctyp) ) 
    stop("invalid ewma type")
  ltyp <- pmatch(limits, c("fix","vacl","fir","both","Steiner","stat","cfar")) - 1
  if ( is.na(ltyp) ) 
    stop("invalid limits type")
  if ( (sided=="one") & !(limits %in% c("fix", "vacl", "stat")) )
    stop("not supported for one-sided EWMA (not reasonable or not implemented yet")
  if ( is.null(c0) ) {
    if ( sided=="one" ) c0 <- zr - 1
    if ( sided=="two" ) c0 <- -1
  }
  if (limits=="cfar") {
    na <- nmax+1
    ctyp <- na
    zero <- function(x) L0 - xewma.ad(l, x, 0, sided="two", r=r)
    c1 <- 1
    if ( zero(c1) > 0 ) {
      while ( zero(c1) > 0 ) c1 <- c1 + 1
      c0 <- c1 - 1
    } else {
      while ( zero(c1) < 0 ) c1 <- c1 / 1.2
      c0 <- c1
      c1 <- c1 * 1.2
    }
    hs <- uniroot(zero, c(c0,c1), tol=1e-9)$root  
    hs <- sqrt( l/(2-l) ) * hs
  } else {
    na <- 1
  }
  cE <- .C("xewma_crit",as.integer(ctyp),as.double(l),
           as.double(L0),as.double(zr),as.double(hs),
           as.double(mu0),as.integer(ltyp),as.integer(r),
           as.double(c0),
           ans=double(length=na),PACKAGE="spc")$ans 
  if (limits=="cfar") {
    nm <- as.integer(cE[1])
    cE <- cE[2:(nm+1)]
    names(cE) <- NULL
  } else {
    names(cE) <- "cE"
  }
  return (cE)
}

