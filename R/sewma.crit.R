# Computation of EWMA critical values for given ARL (variance monitoring)
sewma.crit <- function(l, L0, df, sigma0=1, cl=NULL, cu=NULL, hs=NULL, s2.on=TRUE, sided="upper", mode="fixed", ur=4, r=40, qm=30) {
  
  mitte <- sqrt( 2/df ) * gamma( (df+1)/2 )/ gamma( df/2 )
  if ( is.null(hs) ) {
    if ( s2.on ) { hs <- 1 } else { hs <- mitte }
  }
  
  cu0 <- cl0 <- 0
  if ( l<=0 | l>1 )
    stop("l has to be between 0 and 1")
  if ( L0<1 )
    stop("L0 is too small")
  if ( df<1 )
    stop("df must be positive")
  if ( sigma0<=0 )
    stop("sigma0 must be positive")
  if ( sided=="Rupper" ) {
    if ( is.null(cl) )
      stop("set cl")
    if ( cl<=0 )
      stop("cl must be positive")
    cl0 <- cl
    if ( hs<cl0 )
      stop("hs must not be smaller than cl")
  }
  if ( sided=="Rlower" ) {
    if ( is.null(cu) )
      stop("set cu")
    if ( cu<sigma0 )
      stop(paste("cu must be larger than sigma0 =", sigma0))
    cu0 <- cu
    if ( hs>cu0 )
      stop("hs must not be larger than cu")
  }
  if ( sided=="two" & mode=="fixed" ) {
    if ( is.null(cu) )
      stop("set cu")
    if ( cu<sigma0 )
      stop(paste("cu must be larger than sigma0 =", sigma0))
    cu0 <- cu
    if ( hs>cu0 )
      stop("hs must not be larger than cu")
  }
  s_squared <- as.numeric(s2.on)
  if ( !(s_squared %in% c(0,1)) )
    stop("wrong value for s2.on")
  
  ctyp <- pmatch(sided, c("upper", "Rupper", "two", "Rlower")) - 1
  if (is.na(ctyp))
    stop("invalid ewma type")
  ltyp <- pmatch(mode, c("fixed", "unbiased", "eq.tails", "vanilla")) - 1
  if ( is.na(ltyp) )
    stop("invalid limits type")
  if ( r<10 )
    stop("r is too small")
  if ( qm<10 )
    stop("qm is too small")

  if ( isTRUE(all.equal(1, l)) ) {
    if ( sided=="upper" ) {
      cv <- c(0, qchisq( 1-1/L0, df)/df)
    }
    if ( sided=="Rupper" ) {
      cv <- c(0, Inf)
      warning("not useful for Shewhart type chart")
    }
    if ( sided=="Rlower" ) {
      cv <- c(qchisq( 1/L0, df)/df, Inf)
    }
    if ( sided=="two") {
      if ( mode=="fixed" ) {
        a2 <- 1 - pchisq(df*cu, df)
        if ( a2 > 1/L0 ) {
          cv <- c(0, Inf)
          warning("upper limit too small") 
        } else {
          a1 <- 1/L0 - a2
          cv <- c(qchisq(a1, df)/df, cu)
        }        
      }
      if ( mode=="unbiased") {
        a1 <- 1/(2*L0)
        step <- 1/(2*L0)/10
        one <- 1
        for ( j in 1:8 ) {
          for ( i in 1:11 ) {
            a1 <- a1 + step*one
            k1 <- qchisq(a1, df)
            a2 <- 1/L0 - a1
            k2 <- qchisq(1-a2, df)
            condition <- (k2 - k1)/(log(k2) - log(k1))
            if ( (one*condition) > (one*df) ) break
          }
          step <- step/10
          one <- -one
          if ( abs(condition - df) < 1e-10 ) break
        }
        cv <- c(k1, k2)/df
      }
      if ( mode=="eq.tails" ) {
        a <- 1/(2*L0)
        cv <- qchisq(c(a, 1-a), df)/df
      }
      if ( mode=="vanilla" ) {
        a <- 1/L0
        k2 <- qchisq(1-a, df)
        k1 <- 2*df - k2
        if ( k1 > 0 ) {
          zero <- function(x) ( pchisq( 2*df - qchisq(1-x, df), df) + x ) - 1/L0
          a <- uniroot(zero, c(1/(2*L0), 1/L0), tol=1e-12)$root
          k2 <- qchisq(1-a, df)
          k1 <- 2*df - k2
          cv <- c(k1, k2)/df
        } else {
          cv <- c(0, Inf)
          warning("symmetric design not possible") 
        }
      }      
    }
  } else {
    cv <- .C("sewma_crit",as.integer(ctyp),as.integer(ltyp),as.double(l),
             as.double(L0),as.double(cl0),as.double(cu0),as.double(hs),
             as.double(sigma0),as.integer(df),as.integer(r),as.integer(qm),
             as.double(ur),as.integer(s_squared),
             ans=double(length=2),PACKAGE="spc")$ans
  }
  names(cv) <- c("cl", "cu")
  cv
}

