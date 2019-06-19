\name{tewma.arl}
\alias{tewma.arl}
\title{Compute ARLs of Poisson TEWMA control charts}
\description{Computation of the (zero-state) Average Run Length (ARL) at given Poisson mean \code{mu}.}
\usage{tewma.arl(lambda, k, lk, uk, mu, z0, rando=FALSE, gl=0, gu=0)}
\arguments{
\item{lambda}{smoothing parameter of the EWMA p control chart.}
\item{k}{resolution of grid (natural number).}
\item{lk}{lower control limit of the TEWMA control chart, integer.}
\item{uk}{upper control limit of the TEWMA control chart, integer.}
\item{mu}{mean value of Poisson distribution.}
\item{z0}{so-called headstart (give fast initial response) -- it is proposed to use the in-control mean.}
\item{rando}{Distinguish between control chart design without or with randomisation. In the latter case some
meaningful values for \code{gl} and \code{gu} should be provided.}
\item{gl}{randomisation probability at the lower limit.}
\item{gu}{randomisation probability at the upper limit.}
}
\details{
A new idea of applying EWMA smoothing to count data. Here, the thinning operation is
applied to independent Poisson variates is performed. 
Moreover, the original thinning principle is expanded to multiples of one over \code{k} to allow finer
grids and finally better detection perfomance. It is highly recommended to read the
corresponding paper (see below).
}
\value{Return single value which resemble the ARL.}
\references{
M. C. Morais, C. H. Weiss, S. Knoth (2019),
A thinning-based EWMA chart to monitor counts, submitted.
}
\author{Sven Knoth}
\seealso{later.}
\examples{
# MWK (2018)
lambda <- 0.1 # (T)EWMA smoothing constant
mu0 <- 5 # in-control mean
k <- 10 # resolution
z0 <- round(k*mu0) # starting value of (T)EWMA sequence
# (i) without randomisation
lk <- 28
uk <- 75
L0 <- tewma.arl(lambda, k, lk, uk, mu0, z0)
# should be 501.9703
# (ii) with randomisation
uk <- 76 # lk is not changed
gl <- 0.5446310
gu <- 0.1375617
L0 <- tewma.arl(lambda, k, lk, uk, mu0, z0, rando=TRUE, gl=gl, gu=gu)
# should be 500
}
\keyword{ts}
