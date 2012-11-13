\name{xgrsr.ad}
\alias{xgrsr.ad}
\title{Compute steady-state ARLs of Shiryaev-Roberts schemes}
\description{Computation of the steady-state Average Run Length (ARL)
for Shiryaev-Roberts schemes monitoring normal mean.}
\usage{xgrsr.ad(k, g, mu1, mu0 = 0, zr = 0, sided = "one", MPT = FALSE, r = 30)}
\arguments{
\item{k}{reference value of the Shiryaev-Roberts scheme.}
\item{g}{control limit (alarm threshold) of Shiryaev-Roberts scheme.}
\item{mu1}{out-of-control mean.}
\item{mu0}{in-control mean.}
\item{zr}{reflection border to enable the numerical algorithms used here.}
\item{sided}{distinguishes between one- and two-sided schemes by choosing 
\code{"one"} and\code{"two"}, respectively. Currently only one-sided schemes are
implemented.}
\item{MPT}{switch between the old implementation (\code{FALSE}) and the new one (\code{TRUE}) that considers the completed
likelihood ratio. MPT contains the initials of G. Moustakides, A. Polunchenko and A. Tartakovsky.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1}.}
}
\details{
\code{xgrsr.ad} determines the steady-state Average Run Length (ARL) by numerically
solving the related ARL integral equation by means of the Nystroem method
based on Gauss-Legendre quadrature.
}
\value{Returns a single value which resembles the steady-state ARL.}
\references{
S. Knoth (2006),
The art of evaluating monitoring schemes --
how to measure the performance of control charts?
S. Lenz, H. & Wilrich, P. (ed.),
\emph{Frontiers in Statistical Quality Control 8}, Physica Verlag, Heidelberg, Germany, 74-99.

G. Moustakides, A. Polunchenko, A. Tartakovsky (2009),
Numerical comparison of CUSUM and Shiryaev-Roberts procedures for
detectin changes in distributions,
\emph{Communications in Statistics: Theory and Methods 38}, 3225-3239.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} and \code{xcusum-arl} for zero-state ARL computation of EWMA and CUSUM control charts,
respectively, and \code{xgrsr.arl} for the zero-state ARL.
}
\examples{
## Small study to identify appropriate reflection border to mimic unreflected schemes
k <- .5
g <- log(390)
zrs <- -(0:10)
ZRxgrsr.ad <- Vectorize(xgrsr.ad, "zr")
ads <- ZRxgrsr.ad(k, g, 0, zr=zrs)
data.frame(zrs, ads)

## Table 2 from Knoth (2006)
## original values are
#  mu   arl
#  0    689
#  0.5  30
#  1    8.9
#  1.5  5.1
#  2    3.6
#  2.5  2.8
#  3    2.4
#
k <- .5
g <- log(390)
zr <- -5 # see first example
mus <- (0:6)/2
Mxgrsr.ad <- Vectorize(xgrsr.ad, "mu1")
ads <- round(Mxgrsr.ad(k, g, mus, zr=zr), digits=1)
data.frame(mus, ads)

## Table 4 from Moustakides et al. (2009)
## original values are
# gamma  A        STADD/steady-state ARL
# 50     28.02    4.37
# 100    56.04    5.46
# 500    280.19   8.33
# 1000   560.37   9.64
# 5000   2801.75  12.79
# 10000  5603.7   14.17
Gxgrsr.ad  <- Vectorize("xgrsr.ad", "g")
As <- c(28.02, 56.04, 280.19, 560.37, 2801.75, 5603.7)
gs <- log(As)
theta <- 1
zr <- -6
ads <- round(Gxgrsr.ad(theta/2, gs, theta, zr=zr, r=100), digits=2)
data.frame(As, ads)
}
\keyword{ts}
