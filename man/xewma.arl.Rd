\name{xewma.arl}
\alias{xewma.arl}
\title{Compute ARLs of EWMA control charts}
\description{Computation of the (zero-state) Average Run Length (ARL)
for different types of EWMA control charts monitoring normal mean.}
\usage{xewma.arl(l,c,mu,zr=0,hs=0,sided="one",limits="fix",q=1,r=40)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{mu}{true mean.}
\item{zr}{reflection border for the one-sided chart.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided EWMA control chart 
by choosing \code{"one"} and \code{"two"}, respectively.}
\item{limits}{distinguishes between different control limits behavior.}
\item{q}{change point position. For \eqn{q=1} and
\eqn{\mu=\mu_0} and \eqn{\mu=\mu_1}, the usual
zero-state ARLs for the in-control and out-of-control case, respectively,
are calculated. For \eqn{q>1} and \eqn{\mu!=0} conditional delays, that is,
\eqn{E_q(L-q+1|L\ge q)}, will be determined.
Note that mu0=0 is implicitely fixed.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-sided) or \code{r}
(two-sided).}
}
\details{
In case of the EWMA chart with fixed control limits,
\code{xewma.arl} determines the Average Run Length (ARL) by numerically
solving the related ARL integral equation by means of the Nystroem method
based on Gauss-Legendre quadrature.
If \code{limits} is not \code{"fix"}, then the method presented in Knoth (2003) is utilized.
Note that for one-sided EWMA charts (\code{sided}=\code{"one"}), only
\code{"vacl"} and \code{"stat"} are deployed, while for two-sided ones
(\code{sided}=\code{"two"}) also \code{"fir"}, \code{"both"}
(combination of \code{"fir"} and \code{"vacl"}), and \code{"Steiner"} are implemented.
For details see Knoth (2004).
}
\value{Returns a single value which resembles the ARL.}
\references{
K.-H. Waldmann (1986),
Bounds for the distribution of the run length of geometric moving 
average charts, \emph{Appl. Statist. 35}, 151-158.

S. V. Crowder (1987),
A simple method for studying run-length distributions of
exponentially weighted moving average charts,
\emph{Technometrics 29}, 401-407.

J. M. Lucas and M. S. Saccucci (1990),
Exponentially weighted moving average control schemes: Properties
and enhancements, \emph{Technometrics 32}, 1-12.

S. Chandrasekaran, J. R. English and R. L. Disney (1995), 
Modeling and analysis of EWMA control schemes with variance-adjusted 
control limits, \emph{IIE Transactions 277}, 282-290.

T. R. Rhoads, D. C. Montgomery and C. M. Mastrangelo (1996),
Fast initial response scheme for exponentially weighted moving average
control chart, \emph{Quality Engineering 9}, 317-327.

S. H. Steiner (1999),
EWMA control charts with time-varying control limits and fast initial response, 
\emph{Journal of Quality Technology 31}, 75-86.

S. Knoth (2003),
EWMA schemes with non-homogeneous transition kernels,
\emph{Sequential Analysis 22}, 241-255.

S. Knoth (2004),
Fast initial response features for EWMA Control Charts,
\emph{Statistical Papers 46}, 47-64.
}
\author{Sven Knoth}
\seealso{
\code{xcusum.arl} for zero-state ARL computation of CUSUM control charts
and \code{xewma.ad} for the steady-state ARL.
}
\examples{
## Waldmann (1986), one-sided EWMA
l <- .75
round(xewma.arl(l,2*sqrt((2-l)/l),0,zr=-4*sqrt((2-l)/l)),digits=1)
l <- .5
round(xewma.arl(l,2*sqrt((2-l)/l),0,zr=-4*sqrt((2-l)/l)),digits=1)
## original values are 209.3 and 3907.5 (in Table 2).

## Waldmann (1986), two-sided EWMA with fixed control limits
l <- .75
round(xewma.arl(l,2*sqrt((2-l)/l),0,sided="two"),digits=1)
l <- .5
round(xewma.arl(l,2*sqrt((2-l)/l),0,sided="two"),digits=1)
## original values are 104.0 and 1952 (in Table 1).

## Crowder (1987), two-sided EWMA with fixed control limits
l1 <- .5
l2 <- .05
c <- 2
mu <- (0:16)/4
arl1 <- sapply(mu,l=l1,c=c,sided="two",xewma.arl)
arl2 <- sapply(mu,l=l2,c=c,sided="two",xewma.arl)
round(cbind(mu,arl1,arl2),digits=2)

## original results are (in Table 1)
## 0.00 26.45 127.53
## 0.25 20.12  43.94
## 0.50 11.89  18.97
## 0.75  7.29  11.64
## 1.00  4.91   8.38
## 1.25  3.95*  6.56
## 1.50  2.80   5.41
## 1.75  2.29   4.62
## 2.00  1.94   4.04
## 2.25  1.70   3.61
## 2.50  1.51   3.26
## 2.75  1.37   2.99
## 3.00  1.26   2.76
## 3.25  1.18   2.56
## 3.50  1.12   2.39
## 3.75  1.08   2.26
## 4.00  1.05   2.15  (* -- in Crowder (1987) typo!?)

## Lucas/Saccucci (1990)
## two-sided EWMA

## with fixed limits
l1 <- .5
l2 <- .03
c1 <- 3.071
c2 <- 2.437
mu <- c(0,.25,.5,.75,1,1.5,2,2.5,3,3.5,4,5)
arl1 <- sapply(mu,l=l1,c=c1,sided="two",xewma.arl)
arl2 <- sapply(mu,l=l2,c=c2,sided="two",xewma.arl)
round(cbind(mu,arl1,arl2),digits=2)

## original results are (in Table 3)
## 0.00 500.   500.
## 0.25 255.    76.7
## 0.50  88.8   29.3
## 0.75  35.9   17.6
## 1.00  17.5   12.6
## 1.50   6.53   8.07
## 2.00   3.63   5.99
## 2.50   2.50   4.80
## 3.00   1.93   4.03
## 3.50   1.58   3.49
## 4.00   1.34   3.11
## 5.00   1.07   2.55

## with fir feature
l1 <- .5
l2 <- .03
c1 <- 3.071
c2 <- 2.437
hs1 <- c1/2
hs2 <- c2/2
mu <- c(0,.5,1,2,3,5)
arl1 <- sapply(mu,l=l1,c=c1,hs=hs1,sided="two",limits="fir",xewma.arl)
arl2 <- sapply(mu,l=l2,c=c2,hs=hs2,sided="two",limits="fir",xewma.arl)
round(cbind(mu,arl1,arl2),digits=2)

## original results are (in Table 5)
## 0.0 487.   406.
## 0.5  86.1   18.4
## 1.0  15.9    7.36
## 2.0   2.87   3.43
## 3.0   1.45   2.34
## 5.0   1.01   1.57

## Chandrasekaran, English, Disney (1995)
## two-sided EWMA with fixed and variance adjusted limits (vacl)

l1 <- .25
l2 <- .1
c1s <- 2.9985
c1n <- 3.0042
c2s <- 2.8159
c2n <- 2.8452
mu <- c(0,.25,.5,.75,1,2)
arl1s <- sapply(mu,l=l1,c=c1s,sided="two",xewma.arl)
arl1n <- sapply(mu,l=l1,c=c1n,sided="two",limits="vacl",xewma.arl)
arl2s <- sapply(mu,l=l2,c=c2s,sided="two",xewma.arl)
arl2n <- sapply(mu,l=l2,c=c2n,sided="two",limits="vacl",xewma.arl)
round(cbind(mu,arl1s,arl1n,arl2s,arl2n),digits=2)

## original results are (in Table 2)
## 0.00 500.   500.   500.   500.
## 0.25 170.09 167.54 105.90  96.6
## 0.50  48.14  45.65  31.08  24.35
## 0.75  20.02  19.72  15.71  10.74
## 1.00  11.07   9.37  10.23   6.35
## 2.00   3.59   2.64   4.32   2.73

## The results in Chandrasekaran, English, Disney (1995) are not
## that accurate. Let us consider the more appropriate comparison

c1s <- xewma.crit(l1,500,sided="two")
c1n <- xewma.crit(l1,500,sided="two",limits="vacl")
c2s <- xewma.crit(l2,500,sided="two")
c2n <- xewma.crit(l2,500,sided="two",limits="vacl")
mu <- c(0,.25,.5,.75,1,2)
arl1s <- sapply(mu,l=l1,c=c1s,sided="two",xewma.arl)
arl1n <- sapply(mu,l=l1,c=c1n,sided="two",limits="vacl",xewma.arl)
arl2s <- sapply(mu,l=l2,c=c2s,sided="two",xewma.arl)
arl2n <- sapply(mu,l=l2,c=c2n,sided="two",limits="vacl",xewma.arl)
round(cbind(mu,arl1s,arl1n,arl2s,arl2n),digits=2)

## which demonstrate the abilities of the variance-adjusted limits
## scheme more explicitely.

## Rhoads, Montgomery, Mastrangelo (1996)
## two-sided EWMA with fixed and variance adjusted limits (vacl),
## with fir and both features

l <- .03
c <- 2.437
mu <- c(0,.5,1,1.5,2,3,4)
sl <- sqrt(l*(2-l))
arlfix  <- sapply(mu,l=l,c=c,sided="two",xewma.arl)
arlvacl <- sapply(mu,l=l,c=c,sided="two",limits="vacl",xewma.arl)
arlfir  <- sapply(mu,l=l,c=c,hs=c/2,sided="two",limits="fir",xewma.arl)
arlboth <- sapply(mu,l=l,c=c,hs=c/2*sl,sided="two",limits="both",xewma.arl)
round(cbind(mu,arlfix,arlvacl,arlfir,arlboth),digits=1)

## original results are (in Table 1)
## 0.0 477.3* 427.9* 383.4* 286.2*
## 0.5  29.7   20.0   18.6   12.8
## 1.0  12.5    6.5    7.4    3.6
## 1.5   8.1    3.3    4.6    1.9
## 2.0   6.0    2.2    3.4    1.4
## 3.0   4.0    1.3    2.4    1.0
## 4.0   3.1    1.1    1.9    1.0
## * -- the in-control values differ sustainably from the true values!

## Steiner (1999)
## two-sided EWMA control charts with various modifications

## fixed vs. variance adjusted limits

l <- .05
c <- 3
mu <- c(0,.25,.5,.75,1,1.5,2,2.5,3,3.5,4)
arlfix <- sapply(mu,l=l,c=c,sided="two",xewma.arl)
arlvacl <- sapply(mu,l=l,c=c,sided="two",limits="vacl",xewma.arl)
round(cbind(mu,arlfix,arlvacl),digits=1)

## original results are (in Table 2)
## 0.00 1379.0   1353.0
## 0.25  135.0    127.0
## 0.50   37.4     32.5 
## 0.75   20.0     15.6
## 1.00   13.5      9.0
## 1.50    8.3      4.5
## 2.00    6.0      2.8
## 2.50    4.8      2.0
## 3.00    4.0      1.6
## 3.50    3.4      1.3
## 4.00    3.0      1.1

## fir, both, and Steiner's modification

l <- .03
cfir <- 2.44
cboth <- 2.54
cstein <- 2.55
hsfir <- cfir/2
hsboth <- cboth/2*sqrt(l*(2-l))
mu <- c(0,.5,1,1.5,2,3,4)
arlfir <- sapply(mu,l=l,c=cfir,hs=hsfir,sided="two",limits="fir",xewma.arl)
arlboth <- sapply(mu,l=l,c=cboth,hs=hsboth,sided="two",limits="both",xewma.arl)
arlstein <- sapply(mu,l=l,c=cstein,sided="two",limits="Steiner",xewma.arl)
round(cbind(mu,arlfir,arlboth,arlstein),digits=1)

## original values are (in Table 5)
## 0.0 383.0   384.0   391.0
## 0.5  18.6    14.9    13.8
## 1.0   7.4     3.9     3.6
## 1.5   4.6     2.0     1.8
## 2.0   3.4     1.4     1.3
## 3.0   2.4     1.1     1.0
## 4.0   1.9     1.0     1.0

## SAS/QC manual 1999
## two-sided EWMA control charts with fixed limits

l <- .25
c <- 3
mu <- 1
print(xewma.arl(l,c,mu,sided="two"),digits=11)

# original value is 11.154267016.

## Some recent examples for one-sided EWMA charts
## with varying limits and in the so-called stationary mode

# 1. varying limits = "vacl"

lambda <- .1
L0 <- 500

## Monte Carlo results (10^9 replicates)
# mu    ARL      s.e.
# 0     500.00   0.0160
# 0.5   21.637   0.0006
# 1     6.7596   0.0001
# 1.5   3.5398   0.0001
# 2     2.3038   0.0000
# 2.5   1.7004   0.0000
# 3     1.3675   0.0000

zr <- -6
r <- 50
c <- xewma.crit(lambda, L0, zr=zr, limits="vacl", r=r)
Mxewma.arl <- Vectorize(xewma.arl, "mu")
mus <- (0:6)/2
arls <- round(Mxewma.arl(lambda, c, mus, zr=zr, limits="vacl", r=r), digits=4)
data.frame(mus, arls)

# 2. stationary mode, i. e. limits = "stat"

## Monte Carlo results (10^9 replicates)
# mu    ARL      s.e.
# 0     500.00   0.0159
# 0.5   22.313   0.0006
# 1     7.2920   0.0001
# 1.5   3.9064   0.0001
# 2     2.5131   0.0000
# 2.5   1.7983   0.0000
# 3     1.4029   0.0000

c <- xewma.crit(lambda, L0, zr=zr, limits="stat", r=r)
arls <- round(Mxewma.arl(lambda, c, mus, zr=zr, limits="stat", r=r), digits=4)
data.frame(mus, arls)
}
\keyword{ts}
