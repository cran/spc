\name{mewma.crit}
\alias{mewma.crit}
\title{Compute alarm threshold of MEWMA control charts}
\description{Computation of the alarm threshold for multivariate exponentially weighted
moving average (MEWMA) charts monitoring multivariate normal mean.}
\usage{mewma.crit(l, L0, p, hs=0, r=20)}
\arguments{
\item{l}{smoothing parameter lambda of the MEWMA control chart.}
\item{L0}{in-control ARL.}
\item{p}{dimension of multivariate normal distribution.}
\item{hs}{so-called headstart (enables fast initial response) -- must be non-negative.}
\item{r}{number of quadrature nodes -- dimension of the resulting linear equation system.}
}
\details{
\code{mewma.crit} determines the alarm threshold of for given in-control ARL \code{L0}
by applying secant rule and using \code{mewma.arl()} with \code{ntype="gl2"}.
}
\value{Returns a single value which resembles the critical value \code{c}.}
\references{
Steven E. Rigdon (1995), An integral equation for the in-control average run length of a multivariate
exponentially weighted moving average control chart, \emph{J. Stat. Comput. Simulation 52}, 351-365.
}
\author{Sven Knoth}
\seealso{\code{mewma.arl} for zero-state ARL computation.}
\examples{
# Rigdon (1995), p. 358, Tab. 1
p <- 4
L0 <- 500
r <- .25
h4 <- mewma.crit(r, L0, p)
h4
## original value is 16.38.
}
\keyword{ts}
