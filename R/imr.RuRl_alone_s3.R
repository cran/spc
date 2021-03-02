imr.RuRl_alone_s3 <- function(L0, N=30, qm=30, M0=12) {
  zero <- Vectorize(function(x) imr.arl(M0, x, 0, 1, vsided="two", Rl=1/x^3, N=N, qm=qm) - L0)                    
  Ru1 <- sqrt(2) * qnorm( 1-1/(4*L0) )
  D1 <- zero(Ru1)  
  if ( D1 < 0 ) {
    while ( D1 < 0 ) {
      Ru2 <- Ru1
      Ru1 <- Ru1 * 1.1
      D1 <- zero(Ru1)
    }
    D1 <- Ru1
    Ru1 <- Ru2
    Ru2 <- D1
  } else {
    while ( D1 > 0 ) {
      Ru2 <- Ru1
      Ru1 <- Ru1 / 1.1
      D1 <- zero(Ru1)
    }
  }
  Ru <- uniroot(zero, c(Ru1, Ru2), tol=1e-9)$root
  Rl <- 1/Ru^3
  c(Rl, Ru)
}
