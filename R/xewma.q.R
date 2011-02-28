# Computation of EWMA quantiles (mean monitoring)
xewma.q <- function(l, c, mu, p, zr=0, hs=0, sided="one", r=40) {
  if ( l <= 0 | l > 1 )
    stop("l (lambda) has to be between 0 and 1")
  if ( c <= 0 )
    stop("critical value c has to be positive")
  if ( p <= 0 | p >= 1)
    stop("quantile level p must be in (0,1)")
  if ( zr > c & sided == "one")
    stop("wrong reflexion border")
  if ( (sided == "two" & abs(hs) > c) | (sided == "one" & ( hs < zr | hs > c )) )
    stop("wrong headstart")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if ( is.na(ctyp) )
    stop("invalid ewma type")
  if ( r < 4 )
    stop("r (dimension of Markov chain) is too small")
  quant <- .C("xewma_q", as.integer(ctyp), as.double(l),
              as.double(c), as.double(p), as.double(zr),
              as.double(hs), as.double(mu), as.integer(r),
              ans=double(length=1),PACKAGE="spc")$ans
  names(quant) <- "q"
  quant
}
