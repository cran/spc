\name{tol.lim.fac}
\alias{tol.lim.fac}
\title{Two-sided tolerance limit factors}
\description{For constructing tolerance intervals, which
cover a given proportion \eqn{p}{p} of a normal distribution with
unknown mean and variance with confidence 
\eqn{1-\alpha}{1-a}, one needs to calculate
the so-called tolerance limit factors \eqn{k}{k}. These values
are computed for a given sample size \eqn{n}{n}.}
\usage{tol.lim.fac(n,p,a,mode="WW",m=30)}
\arguments{
\item{n}{sample size.}
\item{p}{coverage.}
\item{a}{error probability \eqn{\alpha}{a}, resulting interval covers at least proportion \code{p}
with confidence of at least \eqn{1-\alpha}{1-a}.}
\item{mode}{distinguish between Wald/Wolfowitz' approximation method (\code{"WW"}) and the more
accurate approach (\code{"exact"})
based on Gauss-Legendre quadrature.}
\item{m}{number of abscissas for the quadrature (needed only for \code{method="exact"}),
of course, the larger the more accurate.}
}
\details{
\code{tol.lim.fac} determines tolerance limits factors 
\eqn{k}{k}
by means of the fast and simple approximation due to
Wald/Wolfowitz (1946) and of Gauss-Legendre quadrature like Odeh/Owen
(1980), respectively, who used in fact the Simpson Rule. Then, by
\eqn{\bar x \pm k \cdot s}{xbar +- k s} 
one can build the tolerance intervals
which cover at least proportion \eqn{p}{p} of a normal distribution for
given confidence level of 
\eqn{1-\alpha}{1-a}. \eqn{\bar x}{xbar} and \eqn{s}{s} stand
for the sample mean and the sample standard deviation, respectively.}
\value{Returns a single value which resembles the tolerance limit factor.}
\references{
A. Wald, J. Wolfowitz (1946), Tolerance limits for a normal distribution,
\emph{Annals of Mathematical Statistics 17}, 208-215.

R. E. Odeh, D. B. Owen (1980), \emph{Tables for Normal Tolerance Limits},
Sampling Plans, and Screening, Marcel Dekker, New York.
}
\author{Sven Knoth}
\seealso{
\code{qnorm} for the ''asymptotic'' case -- cf. second example.
}
\examples{
n <- 2:10
p <- .95
a <- .05
kWW <- sapply(n,p=p,a=a,tol.lim.fac)
kEX <- sapply(n,p=p,a=a,mode="exact",tol.lim.fac)
print(cbind(n,kWW,kEX),digits=4)
## Odeh/Owen (1980), page 98, in Table 3.4.1
##  n  factor k
##  2  36.519
##  3   9.789
##  4   6.341
##  5   5.077
##  6   4.422
##  7   4.020
##  8   3.746
##  9   3.546
## 10   3.393

## n -> infty
n <- 10^{1:7}
p <- .95
a <- .05
kEX <- round(sapply(n,p=p,a=a,mode="exact",tol.lim.fac),digits=4)
kEXinf <- round(qnorm(1-a/2),digits=4)
print(rbind(cbind(n,kEX),c("infinity",kEXinf)),quote=FALSE)
}
\keyword{ts}
