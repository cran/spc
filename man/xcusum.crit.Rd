\name{xcusum.crit}
\alias{xcusum.crit}
\title{Compute decision intervals of CUSUM control charts}
\description{Computation of the decision intervals (alarm limits)
for different types of CUSUM control charts monitoring normal mean.}
\usage{xcusum.crit(k, L0, mu0 = 0, hs = 0, sided = "one", r = 30)}
\arguments{
\item{k}{reference value of the CUSUM control chart.}
\item{L0}{in-control ARL.}
\item{mu0}{in-control mean.}
\item{hs}{so-called headstart (give fast initial response).}
\item{sided}{distinguish between one-, two-sided and Crosier's modified
two-sided CUSUM schemoosing \code{"one"}, \code{"two"},
and \code{"Crosier"}, respectively.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-, two-sided) or \code{2r+1}
(Crosier).}
}
\details{
\code{xcusum.crit} determines the decision interval (alarm limit)
for given in-control ARL \code{L0}
by applying secant rule and using \code{xcusum.arl()}.
}
\value{Returns a single value which resembles the decision interval
\code{h}.}
%\references{Later...}
\author{Sven Knoth}
\seealso{\code{xcusum.arl} for zero-state ARL computation.}
\examples{
k <- .5
incontrolARL <- c(500,5000,50000)
sapply(incontrolARL,k=k,xcusum.crit,r=10) # accuracy with 10 nodes
sapply(incontrolARL,k=k,xcusum.crit,r=20) # accuracy with 20 nodes
sapply(incontrolARL,k=k,xcusum.crit)      # accuracy with 30 nodes
}
\keyword{ts}
