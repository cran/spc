# Computation of EWMA ARLs (mean monitoring)
xewma.arl <- function(l,c,mu,zr=0,hs=0,sided="one",limits="fix",r=40) { 
  if (l<=0 || l>1) 
    stop("l has to be between 0 and 1")
  if (c<=0) 
    stop("c has to be positive")
  if (zr>c & sided=="one") 
    stop("wrong reflexion border")
  if ( (sided=="two" & abs(hs)>c) | (sided=="one" & (hs<zr | hs>c)) ) 
    stop("wrong headstart")
  if (r<4) 
    stop("r is too small")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if (is.na(ctyp)) 
    stop("invalid ewma type")
  ltyp <- pmatch(limits, c("fix","vacl","fir","both","Steiner","Knoth")) - 1
  if (is.na(ltyp)) 
    stop("invalid limits type")
  arl <- .C("xewma_arl",as.integer(ctyp),as.double(l),
            as.double(c),as.double(zr),as.double(hs),
            as.double(mu),as.integer(ltyp),as.integer(r),
            ans=double(length=1),PACKAGE="spc")$ans 
  names(arl) <- "arl"
  return (arl)
}