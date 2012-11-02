\name{xcusum.crit.L0h}
\alias{xcusum.crit.L0h}
\title{Compute the CUSUM reference value k for given in-control ARL and threshold h}
\description{Computation of the reference value k
for one-sided CUSUM control charts monitoring normal mean, if the in-control ARL L0 and
the alarm threshold h are given.}
\usage{xcusum.crit.L0h(L0, h, hs=0, sided="one", r=30, L0.eps=1e-6, k.eps=1e-8)}
\arguments{
\item{L0}{in-control ARL.}
\item{h}{alarm level of the CUSUM control chart.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one-, two-sided and Crosier's modified
two-sided CUSUM scheme choosing \code{"one"}, \code{"two"},
and \code{"Crosier"}, respectively.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-, two-sided) or \code{2r+1}
(Crosier).}
\item{L0.eps}{error bound for the L0 error.}
\item{k.eps}{bound for the difference of two successive values of k.}
}
\details{
\code{xcusum.crit.L0h} determines the reference value k
for given in-control ARL \code{L0} and alarm level \code{h}
by applying secant rule and using \code{xcusum.arl()}. Note that
not for any combination of \code{L0} and \code{h} a solution exists
-- for given \code{L0} there is a maximal value for \code{h}
to get a valid result \code{k}.
}
\value{Returns a single value which resembles the reference value \code{k}.}
%\references{Later...}
\author{Sven Knoth}
\seealso{\code{xcusum.arl} for zero-state ARL computation.}
\examples{
L0 <- 100
h.max <- xcusum.crit(0, L0, 0)
hs <- (300:1)/100
hs <- hs[hs < h.max]
ks <- NULL
for ( h in hs ) ks <- c(ks, xcusum.crit.L0h(L0, h))  
k.max <- qnorm( 1 - 1/L0 )
plot(hs, ks, type="l", ylim=c(0, max(k.max, ks)), xlab="h", ylab="k")
abline(h=c(0, k.max), col="red")
}
\keyword{ts}
