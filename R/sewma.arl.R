# Computation of EWMA ARLs (variance monitoring)
sewma.arl <- function(l,cl,cu,sigma,df,hs=1,sided="upper",r=40,qm=30) { 
  if (l<=0 || l>1) 
    stop("l has to be between 0 and 1")
  if (cu<=0) 
    stop("cu has to be positive")
  if (cl<0)
    stop("cu has to be non-negative")
  if ( hs<cl | hs>cu ) 
    stop("wrong headstart")
  if (r<10) 
    stop("r is too small")
  ctyp <- pmatch(sided, c("upper", "Rupper","two","lower")) - 1
  if (is.na(ctyp)) 
    stop("invalid ewma type")
  arl <- .C("sewma_arl",as.integer(ctyp),as.double(l),
            as.double(cl),as.double(cu),as.double(hs),
            as.double(sigma),as.integer(df),as.integer(r),as.integer(qm),
            ans=double(length=1),PACKAGE="spc")$ans 
  names(arl) <- "arl"
  return (arl)
}
