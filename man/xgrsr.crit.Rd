\name{xgrsr.crit}
\alias{xgrsr.crit}
\title{Compute alarm thresholds for Shiryaev-Roberts schemes}
\description{Computation of the alarm thresholds (alarm limits)
for Shiryaev-Roberts schemes monitoring normal mean.}
\usage{xgrsr.crit(k, L0, mu0 = 0, zr = 0, hs = NULL, sided = "one", r = 30)}
\arguments{
\item{k}{reference value of the Shiryaev-Roberts scheme.}
\item{L0}{in-control ARL.}
\item{mu0}{in-control mean.}
\item{zr}{reflection border to enable the numerical algorithms used here.}
\item{hs}{so-called headstart (enables fast initial response). If \code{hs=NULL}, then
the classical headstart -Inf is used (corresponds to 0 for the non-log scheme).}
\item{sided}{distinguishes between one- and two-sided schemes by choosing 
\code{"one"} and\code{"two"}, respectively. Currently only one-sided schemes are
implemented.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1}.}
}
\details{
\code{xgrsr.crit} determines the alarm threshold (alarm limit)
for given in-control ARL \code{L0}
by applying secant rule and using \code{xgrsr.arl()}.
}
\value{Returns a single value which resembles the alarm limit
\code{g}.}
\references{
G. Moustakides, A. Polunchenko, A. Tartakovsky (2009),
Numerical comparison of CUSUM and Shiryaev-Roberts procedures for detecting changes in distributions,
\emph{Communications in Statistics: Theory and Methods 38}, 3225-3239.r.
}
\author{Sven Knoth}
\seealso{\code{xgrsr.arl} for zero-state ARL computation.}
\examples{
## Table 4 from Moustakides et al. (2009)
## original values are
# gamma/L0  A/exp(g)
# 50        28.02
# 100       56.04
# 500       280.19
# 1000      560.37
# 5000      2801.75
# 10000     5603.7
theta <- 1
zr <- -6
r <- 100
Lxgrsr.crit  <- Vectorize("xgrsr.crit", "L0")
L0s <- c(50, 100, 500, 1000, 5000, 10000)
gs <- Lxgrsr.crit(theta/2, L0s, zr=zr, r=r)
data.frame(L0s, gs, A=round(exp(gs), digits=2))
}
\keyword{ts}
