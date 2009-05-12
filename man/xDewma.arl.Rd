\name{xDewma.arl}
\alias{xDewma.arl}
\title{Compute ARLs of EWMA control charts under drift}
\description{Computation of the (zero-state and other)
Average Run Length (ARL) under drift
for different types of EWMA control charts monitoring normal mean.}
\usage{xDewma.arl(l, c, delta, zr = 0, hs = 0, sided = "one", limits = "fix",
    mode = "Gan", m = NULL, q = 1, r = 40, with0 = FALSE)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{delta}{true drift parameter.}
\item{zr}{reflection border for the one-sided chart.}
\item{hs}{so-called headstart (give fast initial response).}
\item{sided}{distinguish between one- and two-sided EWMA control chart 
by choosing \code{"one"} and \code{"two"}, respectively.}
\item{limits}{distinguish between different control limits behavior.}
\item{mode}{decide whether Gan's or Knoth's approach is used. Use
\code{"Gan"} and \code{"Knoth"}, respectively.}
\item{m}{parameter used if \code{mode="Gan"}. \code{m} is design
parameter of Gan's approach. If \code{m=NULL}, then \code{m}
will increased until the resulting ARL does not change anymore.}
\item{q}{change point position. For \eqn{q=1} and
\eqn{\mu=\mu_0} and \eqn{\mu=\mu_1}, the usual
zero-state ARLs for the in-control and out-of-control case, respectively,
are calculated. For \eqn{q>1} and \eqn{\mu!=0} conditional delays, that is,
\eqn{E_q(L-q+1|L\geq)}, will be determined.
Note that mu0=0 is implicitely fixed. Deploy large \code{q} to mimic
steady-state. It works only for \code{mode="Knoth"}}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-sided) or \code{r}
(two-sided).}
\item{with0}{define whether the first observation used for the ARL
calculatio follows already 1*delta or 0*delta.
With \code{q} additional flexibility is given.}
}
\details{
Based on Gan (1991) or Knoth (2003), the ARL is calculated for EWMA
control charts under drift. In case of Gan's framework, the usual
ARL function with mu=m*delta is determined and recursively via
m-1, m-2, ... 1 (or 0) the drift ARL determined. The framework
of Knoth allows to calculate ARLs for varying parameters, such as
control limits and distributional parameters.
For details see the cited papers.
}
\value{Returns a single value which resembles the ARL.}
\references{
F. F. Gan (1991),
EWMA control chart under linear drift,
\emph{J. Stat. Comput. Simulation 38}, 181-200.

L. A. Aerne, C. W. Champ and S. E. Rigdon (1991),
Evaluation of control charts under linear trend,
\emph{Commun. Stat., Theory Methods 20}, 3341-3349.

S. Knoth (2003),
EWMA schemes with non-homogeneous transition kernels,
\emph{Sequential Analysis 22}, 241-255.

H. M. Fahmy and E. A. Elsayed (2006),
Detection of linear trends in process mean,
\emph{International Journal of Production Research 44}, 487-504.

C. Zou, Y. Liu and Z. Wang (2009),
Comparisons of control schemes for monitoring
the means of processes subject to drifts,
\emph{to appear in Metrika}.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} and \code{xewma.ad} for zero-state and
steady-state ARL computation of EWMA control charts
for the classical step change model.
}
\examples{
DxDewma.arl <- Vectorize(xDewma.arl, "delta")
## Gan (1991)
## Table 1
## original values are
#  delta   arlE1  arlE2  arlE3
#  0       500    500    500
#  0.0001  482    460    424
#  0.0010  289    231    185
#  0.0020  210    162    129
#  0.0050  126     94.6   77.9
#  0.0100   81.7   61.3   52.7
#  0.0500   27.5   21.8   21.9
#  0.1000   17.0   14.2   15.3
#  1.0000    4.09   4.28   5.25
#  3.0000    2.60   2.90   3.43
#
lambda1 <- 0.676
lambda2 <- 0.242
lambda3 <- 0.047
h1 <- 2.204/sqrt(lambda1/(2-lambda1))
h2 <- 1.111/sqrt(lambda2/(2-lambda2))
h3 <- 0.403/sqrt(lambda3/(2-lambda3))
deltas <- c(.0001, .001, .002, .005, .01, .05, .1, 1, 3)
arlE10 <- round(xewma.arl(lambda1, h1, 0, sided="two"), digits=2)
arlE1 <- c(arlE10, round(DxDewma.arl(lambda1, h1, deltas, sided="two", with0=TRUE),
                         digits=2))
arlE20 <- round(xewma.arl(lambda2, h2, 0, sided="two"), digits=2)
arlE2 <- c(arlE20, round(DxDewma.arl(lambda2, h2, deltas, sided="two", with0=TRUE),
                         digits=2))
arlE30 <- round(xewma.arl(lambda3, h3, 0, sided="two"), digits=2)
arlE3 <- c(arlE30, round(DxDewma.arl(lambda3, h3, deltas, sided="two", with0=TRUE),
                         digits=2))
data.frame(delta=c(0, deltas), arlE1, arlE2, arlE3)

## do the same with more digits for the alarm threshold
L0 <- 500
h1 <- xewma.crit(lambda1, L0, sided="two")
h2 <- xewma.crit(lambda2, L0, sided="two")
h3 <- xewma.crit(lambda3, L0, sided="two")
lambdas <- c(lambda1, lambda2, lambda3)
hs <- c(h1, h2, h3) * sqrt(lambdas/(2-lambdas))
hs
# compare with Gan's values 2.204, 1.111, 0.403
round(hs, digits=3)

arlE10 <- round(xewma.arl(lambda1, h1, 0, sided="two"), digits=2)
arlE1 <- c(arlE10, round(DxDewma.arl(lambda1, h1, deltas, sided="two", with0=TRUE),
                         digits=2))
arlE20 <- round(xewma.arl(lambda2, h2, 0, sided="two"), digits=2)
arlE2 <- c(arlE20, round(DxDewma.arl(lambda2, h2, deltas, sided="two", with0=TRUE),
                         digits=2))
arlE30 <- round(xewma.arl(lambda3, h3, 0, sided="two"), digits=2)
arlE3 <- c(arlE30, round(DxDewma.arl(lambda3, h3, deltas, sided="two", with0=TRUE),
                         digits=2))
data.frame(delta=c(0, deltas), arlE1, arlE2, arlE3)

## Aerne et al. (1991) -- two-sided EWMA
## Table I (continued)
## original numbers are
#     delta  arlE1  arlE2  arlE3
#  0.000000  465.0  465.0  465.0
#  0.005623  77.01  85.93  102.68
#  0.007499  64.61  71.78  85.74
#  0.010000  54.20  59.74  71.22
#  0.013335  45.20  49.58  58.90
#  0.017783  37.76  41.06  48.54
#  0.023714  31.54  33.95  39.87
#  0.031623  26.36  28.06  32.68
#  0.042170  22.06  23.19  26.73
#  0.056234  18.49  19.17  21.84
#  0.074989  15.53  15.87  17.83
#  0.100000  13.07  13.16  14.55
#  0.133352  11.03  10.94  11.88
#  0.177828   9.33   9.12   9.71
#  0.237137   7.91   7.62   7.95
#  0.316228   6.72   6.39   6.52
#  0.421697   5.72   5.38   5.37
#  0.562341   4.88   4.54   4.44
#  0.749894   4.18   3.84   3.68
#  1.000000   3.58   3.27   3.07
#
lambda1 <- .133
lambda2 <- .25
lambda3 <- .5
cE1 <- 2.856
cE2 <- 2.974
cE3 <- 3.049
deltas <- 10^(-(18:0)/8)
arlE10 <- round(xewma.arl(lambda1, cE1, 0, sided="two"), digits=2)
arlE1 <- c(arlE10, round(DxDewma.arl(lambda1, cE1, deltas, sided="two"), digits=2))
arlE20 <- round(xewma.arl(lambda2, cE2, 0, sided="two"), digits=2)
arlE2 <- c(arlE20, round(DxDewma.arl(lambda2, cE2, deltas, sided="two"), digits=2))
arlE30 <- round(xewma.arl(lambda3, cE3, 0, sided="two"), digits=2)
arlE3 <- c(arlE30, round(DxDewma.arl(lambda3, cE3, deltas, sided="two"), digits=2))
data.frame(delta=c(0, round(deltas, digits=6)), arlE1, arlE2, arlE3)


## Fahmy/Elsayed (2006) -- two-sided EWMA
## Table 4 (Monte Carlo results, 10^4 replicates, change point at t=51!)
## original numbers are
#   delta     arl  s.e.
#   0.00  365.749  3.598
#   0.10   12.971  0.029
#   0.25    7.738  0.015
#   0.50    5.312  0.009
#   0.75    4.279  0.007
#   1.00    3.680  0.006
#   1.25    3.271  0.006
#   1.50    2.979  0.005
#   1.75    2.782  0.004
#   2.00    2.598  0.005
#
lambda <- 0.1
cE <- 2.7
deltas <- c(.1, (1:8)/4)
arlE1 <- c(round(xewma.arl(lambda, cE, 0, sided="two"), digits=3),
           round(DxDewma.arl(lambda, cE, deltas, sided="two"), digits=3))
arlE51 <- c(round(xewma.arl(lambda, cE, 0, sided="two", q=51), digits=3),
     round(DxDewma.arl(lambda, cE, deltas, sided="two", mode="Knoth", q=51),
           digits=3))
data.frame(delta=c(0, deltas), arlE1, arlE51)

## additional Monte Carlo results with 10^8 replicates
#   delta   arl.q=1   s.e.    arl.q=51  s.e.
#   0.00    368.910   0.036   361.714   0.038
#   0.10     12.986   0.000    12.781   0.000
#   0.25      7.758   0.000     7.637   0.000
#   0.50      5.318   0.000     5.235   0.000
#   0.75      4.285   0.000     4.218   0.000
#   1.00      3.688   0.000     3.628   0.000
#   1.25      3.274   0.000     3.233   0.000
#   1.50      2.993   0.000     2.942   0.000
#   1.75      2.808   0.000     2.723   0.000
#   2.00      2.616   0.000     2.554   0.000

## Zou et al. (2009) -- one-sided EWMA
## Table 1
## original values are
#  delta   arl1  arl2  arl3
#  0           ~ 1730
#  0.0005  317   377   440
#  0.001   215   253   297
#  0.005   83.6  92.6  106
#  0.01    55.6  58.8  66.1
#  0.05    22.6  21.1  22.0
#  0.1     15.5  13.9  13.8
#  0.5     6.65  5.56  5.09
#  1.0     4.67  3.83  3.43
#  2.0     3.21  2.74  2.32
#  3.0     2.86  2.06  1.98
#  4.0     2.14  2.00  1.83
l1 <- 0.03479
l2 <- 0.11125
l3 <- 0.23052
c1 <- 2.711
c2 <- 3.033
c3 <- 3.161
zr <- -6
r  <- 50
deltas <- c(0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1:4)
arl1 <- c(round(xewma.arl(l1, c1, 0, zr=zr, r=r), digits=2),
          round(DxDewma.arl(l1, c1, deltas, zr=zr, r=r), digits=2))
arl2 <- c(round(xewma.arl(l2, c2, 0, zr=zr), digits=2),
          round(DxDewma.arl(l2, c2, deltas, zr=zr, r=r), digits=2))
arl3 <- c(round(xewma.arl(l3, c3, 0, zr=zr, r=r), digits=2),
          round(DxDewma.arl(l3, c3, deltas, zr=zr, r=r), digits=2))
data.frame(delta=c(0, deltas), arl1, arl2, arl3)
}
\keyword{ts}
