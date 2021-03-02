imr.Ru_Rlgiven <- function(Rl, L0, N=30, qm=30, M0=12) {
  zero <- function(x) imr.arl(M0, x, 0, 1, vsided="two", Rl=Rl, N=N, qm=qm) - L0
  Ru1 <- sqrt(2) * qnorm(1-1/(4*L0))
  z1 <- zero(Ru1)
  Ru <- 0
  if ( z1 > 0 ) {
    while ( z1 > 0 ) {
      Ru2 <- Ru1
      Ru1 <- Ru1 / 1.1
      z1 <- zero(Ru1)
    }
  } else {
    if ( zero(2*M0-1e-10) > 0 ) {
      while ( z1 < 0 ) {
        Ru2 <- Ru1
        Ru1 <- Ru1 * 1.1
        z1 <- zero(Ru1)
      }
      z1 <- Ru1
      Ru1 <- Ru2
      Ru2 <- z1
    } else {
      Ru <- Inf
    }
  }
  if ( is.finite(Ru) ) Ru <- uniroot(zero, c(Ru1, Ru2), tol=1e-9)$root
  Ru
}
