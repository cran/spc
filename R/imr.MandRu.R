imr.MandRu <- function(L0, N=30, qm=30) {
  Xarl <- function(M) 1/(2*pnorm(-M))
  MRarl <- function(Ru, N=20, qm=30) ifelse(Ru<100, imr.arl(10, Ru, 0, 1, N=N, qm=qm), Inf)
  M0 <- qnorm(1-1/(2*L0))
  M1 <- M0 + .1
  delta <- Vectorize(function(x) {
    LM <- Xarl(x)
    Ru <- imr.Ru_Mgiven(x, L0, N=N, qm=qm)
    LR <- MRarl(Ru, N=N, qm=qm)
    LM - LR
  })
  D1 <- delta(M1)
  if ( D1 > 0 ) {
    while ( D1 > 0 ) {
      M2 <- M1
      D2 <- D1
      M1 <- mean(c(M0,M1))
      D1 <- delta(M1)
    }
    M <- uniroot(delta, c(M1, M2), tol=1e-9)$root
  } else {
    while ( D1 < 0 ) {
      M2 <- M1
      D2 <- D1
      M1 <- M1 + .1
      D1 <- delta(M1)
    }
    M <- uniroot(delta, c(M2, M1), tol=1e-9)$root
  }
  Ru <- imr.Ru_Mgiven(M, L0, N=N, qm=qm)
  c(M, Ru)
}
