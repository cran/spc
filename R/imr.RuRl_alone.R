imr.RuRl_alone <- function(L0, N=30, qm=30, M0=12, eps=1e-3) {
  zero <- Vectorize(function(x) {
    Ru1 <- imr.Ru_Rlgiven(x, L0, N=N, qm=qm, M0=M0)
    Lp <- imr.arl(M0, Ru1, 0, 1+eps, vsided="two", Rl=x, N=N, qm=qm)
    Lm <- imr.arl(M0, Ru1, 0, 1-eps, vsided="two", Rl=x, N=N, qm=qm)
    DELTA <- (Lp-Lm)/(2*eps)
  })
  pivot <- imr.Rl_Mgiven(M0, L0, N=N, qm=qm)
  Rl1 <- pivot / 2
  D1 <- zero(Rl1)
  Rl <- 1
  if ( D1 < 0 ) {
    while ( D1 < 0 & Rl1 < pivot/1.1 ) {
      Rl2 <- Rl1
      Rl1 <- Rl1 * 1.1
      D1 <- zero(Rl1)
      D1
    }
    if ( D1 > 0 ) {
      D1 <- Rl1
      Rl1 <- Rl2
      Rl2 <- D1
    } else {
      Rl <- 0
    }
  } else {
    while ( D1 > 0 ) {
      Rl2 <- Rl1
      Rl1 <- Rl1 / 1.1
      D1 <- zero(Rl1)
    }
  }
  if ( Rl > 0 ) {
    Rl <- uniroot(zero, c(Rl1, Rl2), tol=1e-9)$root
    Ru <- imr.Ru_Rlgiven(Rl, L0, N=N, qm=qm, M0=M0) 
  } else {
    Ru <- 3*M0 # like infinity
    Rl <- pivot
  }
  c(Rl, Ru)
}
