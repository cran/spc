\name{xcusum.arl}
\alias{xcusum.arl}
\title{Compute ARLs of CUSUM control charts}
\description{Computation of the (zero-state) Average Run Length (ARL)
for different types of CUSUM control charts monitoring normal mean.}
\usage{xcusum.arl(k, h, mu, hs = 0, sided = "one", method = "igl", q = 1, r = 30)}
\arguments{
\item{k}{reference value of the CUSUM control chart.}
\item{h}{decision interval (alarm limit, threshold) of the CUSUM control chart.}
\item{mu}{true mean.}
\item{hs}{so-called headstart (give fast initial response).}
\item{sided}{distinguish between one-, two-sided and Crosier's modified
two-sided CUSUM scheme by choosing \code{"one"}, \code{"two"},
and \code{"Crosier"}, respectively.}
\item{method}{deploy the integral equation (\code{"igl"}) or Markov chain approximation
(\code{"mc"}) method to calculate the ARL (currently only for two-sided CUSUM implemented).}
\item{q}{change point position. For \eqn{q=1} and
\eqn{\mu=\mu_0} and \eqn{\mu=\mu_1}, the usual
zero-state ARLs for the in-control and out-of-control case, respectively,
are calculated. For \eqn{q>1} and \eqn{\mu!=0} conditional delays, that is,
\eqn{E_q(L-q+1|L\ge q)}, will be determined.
Note that mu0=0 is implicitely fixed.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-, two-sided) or \code{2r+1}
(Crosier).}
}
\details{
\code{xcusum.arl} determines the Average Run Length (ARL) by numerically
solving the related ARL integral equation by means of the Nystroem method
based on Gauss-Legendre quadrature.
}
\value{Returns a vector of length \code{q} which resembles the ARL and the sequence of conditional expected delays for
\code{q}=1 and \code{q}>1, respectively.}
\references{
A. L. Goel, S. M. Wu (1971),
Determination of A.R.L. and a contour nomogram for CUSUM charts to
control normal mean, \emph{Technometrics 13}, 221-230.

D. Brook, D. A. Evans (1972),
An approach to the probability distribution of cusum run length,
\emph{Biometrika 59}, 539-548.

J. M. Lucas, R. B. Crosier (1982),
Fast initial response for cusum quality-control schemes: 
Give your cusum a headstart, \emph{Technometrics 24}, 199-205.

L. C. Vance (1986), 
Average run lengths of cumulative sum control charts for controlling 
normal means, \emph{Journal of Quality Technology 18}, 189-193.

K.-H. Waldmann (1986),
Bounds for the distribution of the run length of one-sided and
two-sided CUSUM quality control schemes,
\emph{Technometrics 28}, 61-67.

R. B. Crosier (1986),
A new two-sided cumulative quality control scheme,
\emph{Technometrics 28}, 187-194.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} for zero-state ARL computation of EWMA control charts
and \code{xcusum.ad} for the steady-state ARL.
}
\examples{
## Brook/Evans (1972), one-sided CUSUM
## Their results are based on the less accurate Markov chain approach.

k <- .5
h <- 3
round(c( xcusum.arl(k,h,0), xcusum.arl(k,h,1.5) ),digits=2)

## results in the original paper are L0 = 117.59, L1 = 3.75 (in Subsection 4.3).

## Lucas, Crosier (1982)
## (one- and) two-sided CUSUM with possible headstarts

k <- .5
h <- 4
mu <- c(0,.25,.5,.75,1,1.5,2,2.5,3,4,5)
arl1 <- sapply(mu,k=k,h=h,sided="two",xcusum.arl)
arl2 <- sapply(mu,k=k,h=h,hs=h/2,sided="two",xcusum.arl)
round(cbind(mu,arl1,arl2),digits=2)

## results in the original paper are (in Table 1)
## 0.00 168.   149.
## 0.25  74.2   62.7
## 0.50  26.6   20.1
## 0.75  13.3    8.97
## 1.00   8.38   5.29
## 1.50   4.75   2.86
## 2.00   3.34   2.01
## 2.50   2.62   1.59
## 3.00   2.19   1.32
## 4.00   1.71   1.07
## 5.00   1.31   1.01

## Vance (1986), one-sided CUSUM
## The first paper on using Nystroem method and Gauss-Legendre quadrature
## for solving the ARL integral equation (see as well Goel/Wu, 1971)

k <- 0
h <- 10
mu <- c(-.25,-.125,0,.125,.25,.5,.75,1)
round(cbind(mu,sapply(mu,k=k,h=h,xcusum.arl)),digits=2)

## results in the original paper are (in Table 1 incl. Goel/Wu (1971) results)
##  -0.25  2071.51
##  -0.125  400.28
##   0.0    124.66
##   0.125   59.30
##   0.25    36.71
##   0.50    20.37
##   0.75    14.06
##   1.00    10.75

## Waldmann (1986),
## one- and two-sided CUSUM

## one-sided case

k <- .5
h <- 3
mu <- c(-.5,0,.5)
round(sapply(mu,k=k,h=h,xcusum.arl),digits=2)

## results in the original paper are 1963, 117.4, and 17.35, resp.
## (in Tables 3, 1, and 5, resp.).

## two-sided case

k <- .6
h <- 3
round(xcusum.arl(k,h,-.2,sided="two"),digits=1)  # fits to Waldmann's setup

## result in the original paper is 65.4 (in Table 6).

## Crosier (1986), Crosier's modified two-sided CUSUM
## He introduced the modification and evaluated it by means of
## Markov chain approximation

k <- .5
h <- 3.73
mu <- c(0,.25,.5,.75,1,1.5,2,2.5,3,4,5)
round(cbind(mu,sapply(mu,k=k,h=h,sided="Crosier",xcusum.arl)),digits=2)

## results in the original paper are (in Table 3)
## 0.00 168.
## 0.25  70.7
## 0.50  25.1
## 0.75  12.5
## 1.00   7.92
## 1.50   4.49
## 2.00   3.17
## 2.50   2.49
## 3.00   2.09
## 4.00   1.60
## 5.00   1.22

## SAS/QC manual 1999
## one- and two-sided CUSUM schemes

## one-sided

k <- .25
h <- 8
mu <- 2.5
print(xcusum.arl(k,h,mu),digits=12)
print(xcusum.arl(k,h,mu,hs=.1),digits=12)

## original results are 4.1500836225 and 4.1061588131.

## two-sided

print(xcusum.arl(k,h,mu,sided="two"),digits=12)

## original result is 4.1500826715.
}
\keyword{ts}
