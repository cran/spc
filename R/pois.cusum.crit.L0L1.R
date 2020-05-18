# Computation of Poisson CUSUM alarm threshold, randomization constant and target ooc mean (Ewan/Kemp approach)
pois.cusum.crit.L0L1 <- function(mu0, L0, L1, sided="upper", OUTPUT=FALSE) {
  if ( mu0 < 0 )        stop("mu0 has to be positive")
  if ( L1 <= 1 )        stop("L1 has to be > 1") 
  if ( L0 <= L1 )       stop("L0 has to be > L1")
  if ( !(sided %in% c("upper", "lower")) )      stop("sided has to be either 'upper' or 'lower'")
  
  # helper functions to give the out-of-control mean mu1 for given in-control mean mu0 and reference value k
  k_m01 <- function(mu0, mu1) (mu1 - mu0) / (log(mu1) - log(mu0))
  m1_km0 <- function(mu0, k) {
    zero <- function(x) k - k_m01(mu0,x)
    upper <- mu0 + .5
    while ( zero(upper) > 0 ) upper <- upper + 0.5
    mu1 <- uniroot(zero, c(mu0*1.00000001, upper), tol=1e-9)$root
    mu1
  }
  
  L1_1 <- L1 + 1
  
  m1 <- 10
  k1 <- mu0 * m1
  while ( L1_1 > L1 ) {
    k1 <- k1 + 1
    mu1 <- m1_km0(mu0, k1/m1)
    cv1 <- pois.cusum.crit(mu0, k1, L0, m1, sided=sided, rando=TRUE)
    L1_1 <- pois.cusum.arl(mu1, k1, cv1[1], m1, sided=sided, rando=TRUE, gamma=cv1[2])
    if ( OUTPUT ) cat(paste("k1 =", k1, ",\tmu1 = ", mu1, ",\tL1 =", L1_1, "\n"))
  } 
  
  m1 <- 100
  k1 <- 10 * k1
  while ( L1_1 < L1 ) {
    k1 <- k1 - 1
    mu1 <- m1_km0(mu0, k1/m1)
    cv1 <- pois.cusum.crit(mu0, k1, L0, m1, sided=sided, rando=TRUE)
    L1_1 <- pois.cusum.arl(mu1, k1, cv1[1], m1, sided=sided, rando=TRUE, gamma=cv1[2])
    if ( OUTPUT ) cat(paste("k1 =", k1, ",\tmu1 = ", mu1, ",\tL1 =", L1_1, "\n"))
  }
  
  ff <- max(2, ceiling( 3000 / ( cv1[1]/m1 ) ) / m1)
  m1 <- ff * m1
  k1 <- ff * k1

  while ( L1_1 > L1 ) {
    k1 <- k1 + 1
    mu1 <- m1_km0(mu0, k1/m1)
    cv1 <- pois.cusum.crit(mu0, k1, L0, m1, sided=sided, rando=TRUE)
    L1_1 <- pois.cusum.arl(mu1, k1, cv1[1], m1, sided=sided, rando=TRUE, gamma=cv1[2])
    if ( OUTPUT ) cat(paste("k1 =", k1, ",\tmu1 = ", mu1, ",\tL1 =", L1_1, "\n"))
  }
  
  result <- data.frame(m=m1, km=k1, mu1, k=k1/m1, hm=cv1[1], h=cv1[1]/m1, gamma=cv1[2])
  result
}
