imr.arl.C1987b <- function(M, R, mu, N=80) {
  # original FORTRAN code re-coded in R
  H <- 2*M/N
  NN <- N + 1
  B <- rep(-1, NN)
  A <- matrix(NA, nrow=NN, ncol=NN)
  for ( i in 1:NN ) {
    U <- -M + H*(i-1)
    Bi <- min( M, U+R )
    Ai <- max(-M, U-R )
    N1 <- floor( (Ai+M)/H + 1e-9 ) + 1
    N2 <- floor( (Bi+M)/H + 1e-9 ) + 1
    for ( j in 1:NN ) {
      ARG <- Ai + H*(j-N1)
      A[i,j] <- H * dnorm(ARG-mu)
      if ( j  < N1 ) A[i,j] <- 0
      if ( j == N1 ) A[i,j] <- H/2 * dnorm(Ai-mu)
      if ( j == N2 ) A[i,j] <- H/2 * dnorm(Bi-mu)
      if ( j  > N2 ) A[i,j] <- 0
      if ( j == i  ) A[i,j] <- A[i,j] - 1
    }
  }
  X <- solve(A, B)
  ET <- H/2 * ( dnorm(-M-mu)*X[1] + dnorm(M-mu)*X[NN] ) + 1
  for ( i in 2:N ) {
    ARG <- -M + H*(i-1)
    ET <- ET + H*X[i]*dnorm(ARG-mu)
  }
  ET
}
