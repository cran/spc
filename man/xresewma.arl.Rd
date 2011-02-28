\name{x.res.ewma.arl}
\alias{x.res.ewma.arl}
\alias{s.res.ewma.arl}
\alias{xs.res.ewma.arl}
\alias{xs.res.ewma.pms}
\title{Compute ARLs of EWMA residual control charts}
\description{Computation of the (zero-state) Average Run Length (ARL)
for  EWMA residual control charts monitoring normal mean,
variance, or mean and variance simultaneously. Additionally,
the probability of misleading signals (PMS) is calculated.}
\usage{x.res.ewma.arl(l, c, mu, alpha=0, n=5, hs=0, r=40)

s.res.ewma.arl(l, cu, sigma, mu=0, alpha=0, n=5, hs=1, r=40, qm=30)

xs.res.ewma.arl(lx, cx, ls, csu, mu, sigma, alpha=0, n=5, hsx=0, rx=40, hss=1, rs=40, qm=30)

xs.res.ewma.pms(lx, cx, ls, csu, mu, sigma, type="3", alpha=0, n=5, hsx=0, rx=40, hss=1, rs=40, qm=30)
}
\arguments{
\item{l, lx, ls}{smoothing parameter(s) lambda of the EWMA control chart.}
\item{c, cu, cx, csu}{critical value (similar to alarm limit) of the EWMA control charts.}
\item{mu}{true mean.}
\item{sigma}{true standard deviation.}
\item{alpha}{The AR(1) coefficient -- first order autocorrelation of the original data.}
\item{n}{Batch size.}
\item{hs, hsx, hss}{so-called headstart (give fast initial response).}
\item{r, rx, rs}{number of quadrature nodes or size of collocation base, dimension of the resulting linear
equation system is equal to \code{r} (two-sided).}
\item{qm}{number of nodes for collocation quadratures.}
\item{type}{PMS type, for \code{PMS}="3" (the default) the probability of getting a mean signal despite the variance
changed, and for \code{PMS}="4" the opposite case is dealt with.}
}
\details{
The above list of functions provides the application of algorithms developed for iid data to
the residual case. To be more precise, the underlying model is a sequence of normally
distributed batches with size \code{n} with autocorrelation within the batch and independence between the batches
(see also the references below). It is restricted to the classical EWMA chart types, that
is two-sided for the mean, upper charts for the variance, and all equipped with fixed limits.
The autocorrelation is modeled by an AR(1) process with parameter \code{alpha}. Additionally,
with \code{xs.res.ewma.pms} the probability of misleading signals (PMS) of \code{type} is
calculated. This is offered exclusively in this small collection so that for iid data
this function has to be used too (with \code{alpha=0}).
}
\value{Return single values which resemble the ARL and the PMS, respectively.}
\references{
S. Knoth, M. C. Morais, A. Pacheco, W. Schmid (2009),
Misleading Signals in Simultaneous Residual Schemes for the Mean and Variance of a Stationary Process,
\emph{Commun. Stat., Theory Methods 38}, 2923-2943.

S. Knoth, W. Schmid, A. Schoene (2001),
Simultaneous Shewhart-Type Charts for the Mean and the Variance of a Time Series,
\emph{Frontiers of Statistical Quality Control 6,
A. Lenz, H.-J. & Wilrich, P.-T. (Eds.)}, 6, 61-79.

S. Knoth, W. Schmid (2002)
Monitoring the mean and the variance of a stationary process,
\emph{Statistica Neerlandica 56}, 77-100.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl}, \code{sewma.arl}, and \code{xsewma.arl} as more elaborated functions in the iid case.}
\examples{
## S. Knoth, W. Schmid (2002)

cat("\nFragments of Table 2 (n=5, lambda.1=lambda.2)\n")

lambdas <- c(.5, .25, .1, .05)
L0 <- 500
n <- 5

crit <- NULL
for ( lambda in lambdas ) {
  cs <- xsewma.crit(lambda, lambda, L0, n-1) 
  x.e <- round(cs[1], digits=4)
  names(x.e) <- NULL
  s.e <- round((cs[3]-1) * sqrt((2-lambda)/lambda)*sqrt((n-1)/2), digits=4)
  names(s.e) <- NULL
  crit <- rbind(crit, data.frame(lambda, x.e, s.e))
}


## orinal values are (Markov chain approximation with 50 states)
# lambda x.e    s.e
#   0.50 3.2765 4.6439
#   0.25 3.2168 4.0149
#   0.10 3.0578 3.3376
#   0.05 2.8817 2.9103

print(crit)

cat("\nFragments of Table 4 (n=5, lambda.1=lambda.2=0.1)\n\n")

lambda <- .1
# the algorithm used in Knoth/Schmid is less accurate -- proceed with their values
cx <- x.e <- 3.0578
s.e <- 3.3376
csu <- 1 + s.e * sqrt(lambda/(2-lambda))*sqrt(2/(n-1))

alpha <- .3

a.values <- c((0:6)/4, 2)
d.values <- c(1 + (0:5)/10, 1.75 , 2)

arls <- NULL
for ( delta in d.values ) {
  row <- NULL
  for ( mu in a.values ) {
    arl <- round(xs.res.ewma.arl(lambda, cx, lambda, csu, mu*sqrt(n), delta, alpha=alpha, n=n), digits=2)
    names(arl) <- NULL
    row <- c(row, arl)   
  }
  arls <- rbind(arls, data.frame(t(row)))
}
names(arls) <- a.values
rownames(arls) <- d.values

## orinal values are (now Monte-Carlo with 10^6 replicates)
#          0  0.25   0.5 0.75    1 1.25  1.5    2
#1    502.44 49.50 14.21 7.93 5.53 4.28 3.53 2.65
#1.1   73.19 32.91 13.33 7.82 5.52 4.29 3.54 2.66
#1.2   24.42 18.88 11.37 7.44 5.42 4.27 3.54 2.67
#1.3   13.11 11.83  9.09 6.74 5.18 4.17 3.50 2.66
#1.4    8.74  8.31  7.19 5.89 4.81 4.00 3.41 2.64
#1.5    6.50  6.31  5.80 5.08 4.37 3.76 3.28 2.59
#1.75   3.94  3.90  3.78 3.59 3.35 3.09 2.83 2.40
#2      2.85  2.84  2.80 2.73 2.63 2.51 2.39 2.14

print(arls)

## S. Knoth, M. C. Morais, A. Pacheco, W. Schmid (2009)

cat("\nFragments of Table 5 (n=5, lambda=0.1)\n\n")

d.values <- c(1.02, 1 + (1:5)/10, 1.75 , 2)

arl.x <- arl.s <- arl.xs <- PMS.3 <- NULL
for ( delta in d.values ) {
  arl.x  <- c(arl.x,  round(x.res.ewma.arl(lambda, cx/delta, 0, n=n), digits=3))
  arl.s  <- c(arl.s,  round(s.res.ewma.arl(lambda, csu, delta, n=n), digits=3))
  arl.xs <- c(arl.xs, round(xs.res.ewma.arl(lambda, cx, lambda, csu, 0, delta, n=n), digits=3))
  PMS.3  <- c(PMS.3,  round(xs.res.ewma.pms(lambda, cx, lambda, csu, 0, delta, n=n), digits=6))
}

## orinal values are (Markov chain approximation)
# delta   arl.x   arl.s  arl.xs PMS.3
#  1.02 833.086 518.935 323.324 0.381118
#  1.10 454.101  84.208  73.029 0.145005
#  1.20 250.665  25.871  24.432 0.071024
#  1.30 157.343  13.567  13.125 0.047193
#  1.40 108.112   8.941   8.734 0.035945
#  1.50  79.308   6.614   6.493 0.029499
#  1.75  44.128   3.995   3.942 0.021579
#  2.00  28.974   2.887   2.853 0.018220

cbind(delta=d.values, arl.x, arl.s, arl.xs, PMS.3)

cat("\nFragments of Table 6 (n=5, lambda=0.1)\n\n")

alphas <- c(-0.9, -0.5, -0.3, 0, 0.3, 0.5, 0.9)
deltas <- c(0.05, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2)

PMS.4 <- NULL
for ( ir in 1:length(deltas) ) {
  mu <- deltas[ir]*sqrt(n)
  pms <- NULL
  for ( alpha in alphas ) {
    pms <- c(pms, round(xs.res.ewma.pms(lambda, cx, lambda, csu, mu, 1, type="4", alpha=alpha, n=n), digits=6))
  }
  PMS.4 <- rbind(PMS.4, data.frame(delta=deltas[ir], t(pms)))
}
names(PMS.4) <- c("delta", alphas)
rownames(PMS.4) <- NULL

## orinal values are (Markov chain approximation)
#  delta     -0.9     -0.5     -0.3        0      0.3      0.5      0.9
#   0.05 0.055789 0.224521 0.279842 0.342805 0.391299 0.418915 0.471386
#   0.25 0.003566 0.009522 0.014580 0.025786 0.044892 0.066584 0.192023
#   0.50 0.002994 0.001816 0.002596 0.004774 0.009259 0.015303 0.072945
#   0.75 0.006967 0.000703 0.000837 0.001529 0.003400 0.006424 0.046602
#   1.00 0.005098 0.000402 0.000370 0.000625 0.001589 0.003490 0.039978
#   1.25 0.000084 0.000266 0.000202 0.000300 0.000867 0.002220 0.039773
#   1.50 0.000000 0.000256 0.000120 0.000163 0.000531 0.001584 0.042734
#   2.00 0.000000 0.000311 0.000091 0.000056 0.000259 0.001029 0.054543

PMS.4
}
\keyword{ts}
