
imr.Ru_Mgiven <- function(M, L0, N=30, qm=30) {
  M0 <- qnorm( 1-1/(2*L0) )
  if ( M > M0 ) {
    zero <- function(x) imr.arl(M, x, 0, 1, N=N, qm=qm) - L0
    Ru <- uniroot(zero, c(M0, 10), tol=1e-9)$root
  } else {
    Ru <- Inf
  }
  Ru
}
