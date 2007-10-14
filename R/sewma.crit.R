# Computation of EWMA critical values for given ARL (variance monitoring)
sewma.crit <- function(l,L0,sigma0=1,cu=NULL,hs=1,df,s2.on=TRUE,
                       sided="upper",mode="fixed",r=40,qm=30) {
  if (l<=0 || l>1) 
    stop("l has to be between 0 and 1")
  if (L0<1) 
    stop("L0 is too small")
  if (sigma0<=0)
    stop("sigma0 must be positive")
  if (mode=="fixed" & sided=="two") {
    if (is.null(cu)) stop("set cu")
    if (cu<sigma0) stop("cu is too small")
    if (cu<=0) stop("cu must be positive")
    if (hs>cu) stop("hs must be smaller than cu")
    cu0 <- cu
  } else {
    cu0 <- 0
  }
  if (df<1)
    stop("df must be positive")
  s_squared <- as.numeric(s2.on)
  if ( !(s_squared %in% c(0,1)) )
    stop("wrong value for s2.on")
  ctyp <- pmatch(sided, c("upper","Rupper","two","lower")) - 1
  if (is.na(ctyp)) 
    stop("invalid ewma type")
  ltyp <- pmatch(mode, c("fixed","unbiased")) - 1
  if (is.na(ltyp)) 
    stop("invalid limits type")
  if (r<10) 
    stop("r is too small")
  if (qm<10) 
    stop("qm is too small")
  c <- .C("sewma_crit",as.integer(ctyp),as.integer(ltyp),as.double(l),
          as.double(L0),as.double(cu0),as.double(hs),
          as.double(sigma0),as.integer(df),as.integer(r),as.integer(qm),
          as.integer(s_squared),
          ans=double(length=2),PACKAGE="spc")$ans 
  names(c) <- c("cl","cu")
  return (c)
}

