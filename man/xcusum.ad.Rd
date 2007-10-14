\name{xcusum.ad}
\alias{xcusum.ad}
\title{Compute steady-state ARLs of CUSUM control charts}
\description{Computation of the steady-state Average Run Length (ARL)
for different types of CUSUM control charts monitoring normal mean.}
\usage{xcusum.ad(k, h, mu1, mu0 = 0, sided = "one", r = 30)}
\arguments{
\item{k}{reference value of the CUSUM control chart.}
\item{h}{decision interval (alarm limit, threshold) of the CUSUM control chart.}
\item{mu1}{out-of-control mean.}
\item{mu0}{in-control mean.}
\item{sided}{distinguish between one-, two-sided and Crosier's modified
two-sided CUSUM scheme by choosing \code{"one"}, \code{"two"},
and \code{"Crosier"}, respectively.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-, two-sided) or \code{2r+1}
(Crosier).}
}
\details{
\code{xcusum.ad} determines the steady-state Average Run Length (ARL)
by numerically solving the related ARL integral equation by means
of the Nystroem method based on Gauss-Legendre quadrature
and using the power method for deriving the largest in magnitude
eigenvalue and the related left eigenfunction.
}
\value{Returns a single value which resembles the steady-state ARL.}
\references{
R. B. Crosier (1986),
A new two-sided cumulative quality control scheme,
\emph{Technometrics 28}, 187-194.
}
\note{Be cautious in increasing the dimension parameter \code{r} for
two-sided CUSUM schemes. The resulting matrix dimension is \code{r^2} times
\code{r^2}. Thus, go beyond 30 only on fast machines. This is the only case,
were the package routines are based on the Markov chain approach. Moreover,
the two-sided CUSUM scheme needs a two-dimensional Markov chain.} 
\author{Sven Knoth}
\seealso{
\code{xcusum.arl} for zero-state ARL computation and
\code{xewma.ad} for the steady-state ARL of EWMA control charts.
}
\examples{
## comparison of zero-state (= worst case ) and steady-state performance
## for one-sided CUSUM control charts

k <- .5
h <- xcusum.crit(k,500)
mu <- c(0,.5,1,1.5,2)
arl <- sapply(mu,k=k,h=h,xcusum.arl)
ad <- sapply(mu,k=k,h=h,xcusum.ad)
round(cbind(mu,arl,ad),digits=2)

## Crosier (1986), Crosier's modified two-sided CUSUM
## He introduced the modification and evaluated it by means of
## Markov chain approximation

k <- .5
h2 <- 4
hC <- 3.73
mu <- c(0,.25,.5,.75,1,1.5,2,2.5,3,4,5)
ad2 <- sapply(mu,k=k,h=h2,sided="two",r=20,xcusum.ad)
adC <- sapply(mu,k=k,h=hC,sided="Crosier",xcusum.ad)
round(cbind(mu,ad2,adC),digits=2)

## results in the original paper are (in Table 5)
## 0.00 163.   164.
## 0.25  71.6   69.0
## 0.50  25.2   24.3
## 0.75  12.3   12.1
## 1.00   7.68   7.69
## 1.50   4.31   4.39
## 2.00   3.03   3.12
## 2.50   2.38   2.46
## 3.00   2.00   2.07
## 4.00   1.55   1.60
## 5.00   1.22   1.29.
}
\keyword{ts}
