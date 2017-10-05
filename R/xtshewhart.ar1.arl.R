# Computation of the ARL for modified Shewhart charts, AR(1) data
xtshewhart.ar1.arl <- function(alpha, cS, df, delta=0, N1=50, N2=30, N3=2*N2, INFI=10, mode="tan") {
  if ( abs(alpha) >= 1 )    stop("alpha has to be between -1 and 1")
  if ( cS <= 0 )            stop("cS has to be positive")
  ntyp <- -1 + pmatch(mode, c("identity", "sin", "sinh", "tan"))
  if ( is.na(ntyp) )
    stop("substitution type not provided (yet) or simply wrong") 
  arl <- .C("tshewhart_ar1_arl",
              as.double(alpha), as.double(cS), as.double(delta), as.integer(df),
              as.integer(N1), as.integer(N2), as.integer(N3), as.double(INFI), as.integer(ntyp), ans=double(length=1), PACKAGE="spc")$ans
  names(arl) <- NULL
  arl
}
