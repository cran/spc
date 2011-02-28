\name{xewma.crit}
\alias{xewma.crit}
\title{Compute critical values of EWMA control charts}
\description{Computation of the critical values (similar to alarm limits)
for different types of EWMA control charts monitoring normal mean.}
\usage{xewma.crit(l,L0,mu0=0,zr=0,hs=0,sided="one",limits="fix",r=40,c0=NULL)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{L0}{in-control ARL.}
\item{mu0}{in-control mean.}
\item{zr}{reflection border for the one-sided chart.}
\item{hs}{so-called headstart (give fast initial response).}
\item{sided}{distinguish between one- and two-sided
two-sided EWMA control chart by choosing \code{"one"} and \code{"two"},
respectively.}
\item{limits}{distinguish between different control limits behavior.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-sided) or \code{r}
(two-sided).}
\item{c0}{starting value for iteration rule.}
}
\details{
\code{xewma.crit} determines the critical values (similar to alarm limits)
for given in-control ARL \code{L0}
by applying secant rule and using \code{xewma.arl()}.
}
\value{Returns a single value which resembles the critical value
\code{c}.}
\references{
S. V. Crowder (1989), Design of exponentially weighted moving average
schemes, \emph{Journal of Quality Technology 21}, 155-162.
}
\author{Sven Knoth}
\seealso{\code{xewma.arl} for zero-state ARL computation.}
\examples{
l <- .1
incontrolARL <- c(500,5000,50000)
sapply(incontrolARL,l=l,sided="two",xewma.crit,r=35) # accuracy with 35 nodes
sapply(incontrolARL,l=l,sided="two",xewma.crit)      # accuracy with 40 nodes
sapply(incontrolARL,l=l,sided="two",xewma.crit,r=50) # accuracy with 50 nodes

## Crowder (1989)
## two-sided EWMA control charts with fixed limits

l <- c(.05,.1,.15,.2,.25)
L0 <- 250
round(sapply(l,L0=L0,sided="two",xewma.crit),digits=2)

## original values are 2.32, 2.55, 2.65, 2.72, and 2.76.
}
\keyword{ts}
