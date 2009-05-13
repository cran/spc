\name{sewma.arl}
\alias{sewma.arl}
\title{Compute ARLs of EWMA control charts (variance charts)}
\description{Computation of the (zero-state) Average Run Length (ARL)
for different types of EWMA control charts 
(based on the sample variance \eqn{S^2})
monitoring normal variance.}
\usage{sewma.arl(l,cl,cu,sigma,df,s2.on=TRUE,hs=1,sided="upper",r=40,qm=30)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{cl}{lower control limit of the EWMA control chart.}
\item{cu}{upper control limit of the EWMA control chart.}
\item{sigma}{true standard deviation.}
\item{df}{actual degrees of freedom, corresponds to batch size.}
\item{s2.on}{distinguish between \eqn{S^2}{S^2} and \eqn{S}{S} chart.}
\item{hs}{so-called headstart (give fast initial response).}
\item{sided}{distinguish between one- and two-sided 
two-sided EWMA-\eqn{S^2}{S^2} control charts 
by choosing \code{"upper"} (upper chart
without reflection at \code{cl} -- the actual value of
\code{cl} is not used),
\code{"Rupper"} (upper chart with reflection at \code{cl}),
\code{"Rlower"} (lower chart with reflection at \code{cu}),
and \code{"two"} (two-sided chart), 
respectively.} 
\item{r}{dimension of the resulting linear equation system.}
\item{qm}{number of quadrature nodes.}
}
\details{
\code{sewma.arl} determines the Average Run Length (ARL) by numerically
solving the related ARL integral equation by means of 
collocation (Chebyshev polynomials).}
\value{Returns a single value which resembles the ARL.}
\references{
S. Knoth (2005),
Accurate ARL computation for EWMA-\eqn{S^2}{S^2} control charts,
\emph{Statistics and Computing 15}, 341-352.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} for zero-state ARL computation of EWMA control charts
for monitoring normal mean.
}
\examples{
## Knoth (2005)
## compare with Table 1 (p. 347): 249.9997
## Monte Carlo with 10^9 replicates: 249.9892 +/- 0.008
l <- .025
df <- 1
cu <- 1 + 1.661865*sqrt(l/(2-l))*sqrt(2/df) 
sewma.arl(l,0,cu,1,df)
}
\keyword{ts}
