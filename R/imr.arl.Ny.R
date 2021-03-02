imr.arl.Ny <- function(M, R, mu, N, kind="gl") {
  phij <- function(x,y) ifelse(abs(y-x)<=R, dnorm(y, mean=mu), 0)
  one <- rep(1, N)
  I <- diag(1, N)
  if ( kind == "gl" ) {
    GQ <- quadrature.nodes.weights(N, x1=-M, x2=M)
    z <- GQ$nodes
    w <- GQ$weights
  }
  if ( kind == "rectangular" ) {
    b <- 2*M/N
    z <- -M + (0:(N-1) + .5)*b
    w <- rep(b, N)
  }
  if ( kind %in% c("trapezoid", "simpson", "simpson3_8") ) {
    b <- 2*M/(N-1)
    z <- -M + (0:(N-1))*b
    w <- 1:N
  }
  if ( kind == "trapezoid" ) {
    w <- rep(2, N)
    w[c(1, N)] <- 1
    w <- w * b/2    
  }
  if ( kind == "simpson" ) {
    w <- 4*((w %% 2) == 0) + 2*((w %% 2) == 1)
    w[c(1, N)] <- 1
    w <- w * b/3    
  }
  if ( kind == "simpson3_8" ) {
    w <- 3*((w %% 3) != 1) + 2*((w %% 3) == 1)
    w[c(1, N)] <- 1
    w <- w * 3*b/8    
  }
  Q <- outer(z,z,phij)*t(array(w,c(N,N)))
  vLu <- solve(I-Q, one)
  Lu <- Vectorize(function(x) { 1 + ( phij(x,z) * w ) %*% vLu })
  integrand <- function(x) Lu(x) * dnorm(x, mean=mu)
  GQ2 <- quadrature.nodes.weights(200, x1=-M, x2=M)
  z2 <- GQ2$nodes
  w2 <- GQ2$weights
  ET <- 1 + sum( w2 * integrand(z2) )
  ET
}
