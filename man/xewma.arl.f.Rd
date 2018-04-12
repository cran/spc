\name{xewma.arl.f}
\alias{xewma.arl.f}
\title{Compute ARL function of EWMA control charts}
\description{Computation of the (zero-state) Average Run Length (ARL) function for different types of EWMA control charts monitoring normal mean.}
\usage{xewma.arl.f(l,c,mu,zr=0,sided="one",limits="fix",r=40)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{mu}{true mean.}
\item{zr}{reflection border for the one-sided chart.}
\item{sided}{distinguishes between one- and two-sided EWMA control chart by choosing \code{"one"} and \code{"two"}, respectively.}
\item{limits}{distinguishes between different control limits behavior.}
\item{r}{number of quadrature nodes, dimension of the resulting linear equation system is equal to \code{r+1} (one-sided) or \code{r} (two-sided).}
}
\details{
It is a convenience function to yield the ARL as function of the head start \code{hs}. For more details see \code{xewma.arl}.
}
\value{It returns a function of a single argument, \code{hs=x} which maps the head-start value \code{hs} to the ARL.}
\references{
S. V. Crowder (1987),
A simple method for studying run-length distributions of exponentially weighted moving average charts,
\emph{Technometrics 29}, 401-407.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} for zero-state ARL for one specific head-start \code{hs}.
}
\examples{
# will follow
}
\keyword{ts}
