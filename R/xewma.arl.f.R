# Computation of EWMA ARLs (mean monitoring), returns function
xewma.arl.f <- function(l, c, mu, zr=0, sided="one", limits="fix", r=40) {
  if ( l<=0 | l>1 )             stop("l has to be between 0 and 2")
  if ( c<=0 )                   warning("usually, c has to be positive")
  if ( zr>c & sided=="one" )    stop("wrong reflexion border")
  if ( r<4 )                    stop("r is too small")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if ( is.na(ctyp) )            stop("invalid ewma type")
  ltyp <- -1 + pmatch(limits, c("fix", "vacl", "fir", "both", "Steiner", "stat", "fink", "limit", "fixW", "fixC"))
  if ( is.na(ltyp) )            stop("invalid limits type")
  if ( (sided=="one") & !(limits %in% c("fix", "vacl", "stat", "limit", "fixW")) )
                                stop("not supported for one-sided EWMA (not reasonable or not implemented yet")

  LENGTH <- 3*r
  zeug <- .C("xewma_arl_f", as.integer(ctyp), as.double(l), as.double(c), as.double(zr),
                            as.double(mu), as.integer(ltyp), as.integer(r),
                            ans=double(length=LENGTH), PACKAGE="spc")$ans
                            
  g <- zeug[1:r]                          
  w <- zeug[1:r + r]
  z <- zeug[1:r + 2*r]  
   
  arl <- Vectorize( function(x) 1 + sum( w * dnorm( ( z - (1-l)*x ) / l - mu)/l * g )  )
  
  arl   
}
