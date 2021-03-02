imr.Rl_Mgiven <- function(M, L0, N=30, qm=30) {
  M0 <- qnorm( 1-1/(2*L0) )
  if ( M > M0 ) {
    zero <- function(x) imr.arl(M, 3*M, 0, 1, vsided="two", Rl=x, N=N, qm=qm) - L0
    Rl <- uniroot(zero, c(1e-6, 1), tol=1e-9)$root
  } else {
    Rl <- 0
  }
  Rl
}
