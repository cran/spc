\name{xewma.ad}
\alias{xewma.ad}
\title{Compute steady-state ARLs of EWMA control charts}
\description{Computation of the steady-state Average Run Length (ARL)
for different types of EWMA control charts monitoring normal mean.}
\usage{xewma.ad(l,c,mu1,mu0=0,zr=0,sided="one",limits="fix",r=40)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{mu1}{in-control mean.}
\item{mu0}{out-of-control mean.}
\item{zr}{reflection border for the one-sided chart.}
\item{sided}{distinguish between one- and two-sided
two-sided EWMA control chart by choosing \code{"one"} and \code{"two"},
respectively.}
\item{limits}{distinguish between different control limits behavior.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-sided) or \code{r}
(two-sided).}
}
\details{
\code{xewma.ad} determines the steady-state Average Run Length (ARL)
by numerically solving the related ARL integral equation by means
of the Nyström method based on Gauss-Legendre quadrature
and using the power method for deriving the largest in magnitude
eigenvalue and the related left eigenfunction.
}
\value{Returns a single value which resembles the steady-state ARL.}
\references{
S. V. Crowder (1987),
A simple method for studying run-length distributions of
exponentially weighted moving average charts,
\emph{Technometrics 29}, 401-407.

J. M. Lucas, M. S. Saccucci (1990),
Exponentially weighted moving average control schemes: Properties
and enhancements, \emph{Technometrics 32}, 1-12.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} for zero-state ARL computation and
\code{xcusum.ad} for the steady-state ARL of CUSUM control charts.
}
\examples{
## comparison of zero-state (= worst case ) and steady-state performance
## for two-sided EWMA control charts

l <- .1
c <- xewma.crit(l,500,sided="two")
mu <- c(0,.5,1,1.5,2)
arl <- sapply(mu,l=l,c=c,sided="two",xewma.arl)
ad <- sapply(mu,l=l,c=c,sided="two",xewma.ad)
round(cbind(mu,arl,ad),digits=2)

## Lucas/Saccucci (1990)
## Lucas/Saccucci (1990)
## two-sided EWMA

## with fixed limits
l1 <- .5
l2 <- .03
c1 <- 3.071
c2 <- 2.437
mu <- c(0,.25,.5,.75,1,1.5,2,2.5,3,3.5,4,5)
ad1 <- sapply(mu,l=l1,c=c1,sided="two",xewma.ad)
ad2 <- sapply(mu,l=l2,c=c2,sided="two",xewma.ad)
round(cbind(mu,ad1,ad2),digits=2)

## original results are (in Table 3)
## 0.00 499.   480.  
## 0.25 254.    74.1
## 0.50  88.4   28.6
## 0.75  35.7   17.3
## 1.00  17.3   12.5
## 1.50   6.44   8.00
## 2.00   3.58   5.95
## 2.50   2.47   4.78
## 3.00   1.91   4.02
## 3.50   1.58   3.49
## 4.00   1.36   3.09
## 5.00   1.10   2.55.
}
\keyword{ts}
