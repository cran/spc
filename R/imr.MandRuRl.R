imr.MandRuRl <- function(L0, N=30, qm=30) {
  Xarl <- function(M) 1/(2*pnorm(-M))
  MRarl <- function(Ru, Rl, N=N, qm=qm) imr.arl(12, Ru, 0, 1, vsided="two", Rl=Rl, N=N, qm=qm)
  M0 <- qnorm(1-1/(2*L0))
  M1 <- M0 + .1
  delta <- Vectorize(function(x) {
    LM <- Xarl(x)    
    RRR <- imr.RuRl_alone(L0, N=N, qm=qm, M0=x)    
    LR <- MRarl(RRR[2], RRR[1], N=N, qm=qm)
    LM - LR
  })
  D1 <- delta(M1)
  if ( D1 > 0 ) {
    while ( D1 > 0 ) {
      M2 <- M1
      M1 <- mean(c(M0,M1))
      D1 <- delta(M1)
    }
  } else {
    while ( D1 < 0 ) {
      M2 <- M1
      M1 <- M1 + 0.1
      D1 <- delta(M1)
    }
    while ( D1 > 0 ) {
      M2 <- M1
      M1 <- M1 - 0.02
      D1 <- delta(M1)
    }
  }
  M <- uniroot(delta, c(M1, M2), tol=1e-9)$root
  RRR <- imr.RuRl_alone(L0, N=N, qm=qm, M0=M)
  Rl <- RRR[1]
  Ru <- RRR[2]
  c(M, Rl, Ru)
}
