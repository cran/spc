imr.Rl_Rugiven <- function(Ru, L0, N=30, qm=30, M0=12) {
  zero <- function(x) imr.arl(M0, Ru, 0, 1, vsided="two", Rl=x, N=N, qm=qm) - L0
  Rl1 <- sqrt(2) * qnorm(0.5+1/(4*L0))
  z1 <- zero(Rl1)
  if ( z1 < 0 ) {
    while ( z1 < 0 ) {
      Rl2 <- Rl1
      Rl1 <- Rl1 / 1.1
      z1 <- zero(Rl1)
    }
    z1 <- Rl1
    Rl1 <- Rl2
    Rl2 <- z1
  } else {
    while ( z1 > 0 ) {
      Rl2 <- Rl1
      Rl1 <- Rl1 * 1.1
      z1 <- zero(Rl1)
    }
  }
  Rl <- uniroot(zero, c(Rl1, Rl2), tol=1e-9)$root
  Rl
}
