\name{xgrsr.arl}
\alias{xgrsr.arl}
\title{Compute (zero-state) ARLs of Shiryaev-Roberts schemes}
\description{Computation of the (zero-state) Average Run Length (ARL)
for Shiryaev-Roberts schemes monitoring normal mean.}
\usage{xgrsr.arl(k, g, mu, zr = 0, hs=NULL, sided = "one", q = 1, r = 30)}
\arguments{
\item{k}{reference value of the Shiryaev-Roberts scheme.}
\item{g}{control limit (alarm threshold) of Shiryaev-Roberts scheme.}
\item{mu}{true mean.}
\item{zr}{reflection border to enable the numerical algorithms used here.}
\item{hs}{so-called headstart (enables fast initial response). If \code{hs=NULL}, then
the classical headstart -Inf is used (corresponds to 0 for the non-log scheme).}
\item{sided}{distinguishes between one- and two-sided schemes by choosing 
\code{"one"} and\code{"two"}, respectively. Currently only one-sided schemes are
implemented.}
\item{q}{change point position. For \eqn{q=1} and
\eqn{\mu=\mu_0} and \eqn{\mu=\mu_1}, the usual
zero-state ARLs for the in-control and out-of-control case, respectively,
are calculated. For \eqn{q>1} and \eqn{\mu!=0} conditional delays, that is,
\eqn{E_q(L-q+1|L\ge q)}, will be determined.
Note that mu0=0 is implicitely fixed.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1}.}
}
\details{
\code{xgrsr.arl} determines the Average Run Length (ARL) by numerically
solving the related ARL integral equation by means of the Nystroem method
based on Gauss-Legendre quadrature.
}
\value{Returns a single value which resembles the ARL.}
\references{
S. Knoth (2006),
The art of evaluating monitoring schemes --
how to measure the performance of control charts?
S. Lenz, H. & Wilrich, P. (ed.),
\emph{Frontiers in Statistical Quality Control 8}, Physica Verlag, Heidelberg, Germany, 74-99.

G. Moustakides, A. Polunchenko, A. Tartakovsky (2009),
Numerical comparison of CUSUM and Shiryaev-Roberts procedures for detecting changes in distributions,
\emph{Communications in Statistics: Theory and Methods 38}, 3225-3239.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} and \code{xcusum-arl} for zero-state ARL computation of EWMA and CUSUM control charts,
respectively, and \code{xgrsr.ad} for the steady-state ARL.
}
\examples{
## Small study to identify appropriate reflection border to mimic unreflected schemes
k <- .5
g <- log(390)
zrs <- -(0:10)
ZRxgrsr.arl <- Vectorize(xgrsr.arl, "zr")
arls <- ZRxgrsr.arl(k, g, 0, zr=zrs)
data.frame(zrs, arls)

## Table 2 from Knoth (2006)
## original values are
#  mu   arl
#  0    697
#  0.5  33
#  1    10.4
#  1.5  6.2
#  2    4.4
#  2.5  3.5
#  3    2.9
#
k <- .5
g <- log(390)
zr <- -5 # see first example
mus <- (0:6)/2
Mxgrsr.arl <- Vectorize(xgrsr.arl, "mu")
arls <- round(Mxgrsr.arl(k, g, mus, zr=zr), digits=1)
data.frame(mus, arls)

## Table 4 from Moustakides et al. (2009)
## original values are
# gamma  A        ARL/E_infty(L) SADD/E_1(L)
# 50     28.02    50.79          5.46
# 100    56.04    100.79         6.71
# 500    280.19   500.8          9.78
# 1000   560.37   1000.79        11.14
# 5000   2801.75  5001.75        14.34
# 10000  5603.7   10000.78       15.73
Gxgrsr.arl  <- Vectorize("xgrsr.arl", "g")
As <- c(28.02, 56.04, 280.19, 560.37, 2801.75, 5603.7)
gs <- log(As)
theta <- 1
zr <- -6
arls0 <- round(Gxgrsr.arl(theta/2, gs, 0, zr=zr, r=100), digits=2)
arls1 <- round(Gxgrsr.arl(theta/2, gs, theta, zr=zr, r=100), digits=2)
data.frame(As, arls0, arls1)
}
\keyword{ts}
