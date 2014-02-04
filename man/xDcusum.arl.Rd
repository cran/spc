\name{xDcusum.arl}
\alias{xDcusum.arl}
\title{Compute ARLs of CUSUM control charts under drift}
\description{Computation of the (zero-state and other) Average Run Length (ARL)
under drift for one-sided CUSUM control charts monitoring normal mean.}
\usage{xDcusum.arl(k, h, delta, hs = 0, sided = "one",
    mode = "Gan", m = NULL, q = 1, r = 30, with0 = FALSE)}
\arguments{
\item{k}{reference value of the CUSUM control chart.}
\item{h}{decision interval (alarm limit, threshold) of the CUSUM control chart.}
\item{delta}{true drift parameter.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided CUSUM control chart
by choosing \code{"one"} and \code{"two"}, respectively. Currentlly,
the two-sided scheme is not implemented.}
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
steady-state. It works only for \code{mode="Knoth"}.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-sided) or \code{r}
(two-sided).}
\item{with0}{defines whether the first observation used for the RL
calculation follows already 1*delta or still 0*delta.
With \code{q} additional flexibility is given.}
}
\details{
Based on Gan (1991) or Knoth (2003), the ARL is calculated for
one-sided CUSUM
control charts under drift. In case of Gan's framework, the usual
ARL function with mu=m*delta is determined and recursively via
m-1, m-2, ... 1 (or 0) the drift ARL determined. The framework
of Knoth allows to calculate ARLs for varying parameters, such as
control limits and distributional parameters.
For details see the cited papers. Note that two-sided CUSUM charts
under drift are difficult to treat.
}
\value{Returns a single value which resembles the ARL.}
\references{
F. F. Gan (1992),
CUSUM control charts under linear drift,
\emph{Statistician 41}, 71-84.

F. F. Gan (1996),
Average Run Lengths for Cumulative Sum control chart under linear trend,
\emph{Applied Statistics 45}, 505-512.

S. Knoth (2003),
EWMA schemes with non-homogeneous transition kernels,
\emph{Sequential Analysis 22}, 241-255.

S. Knoth (2012),
More on Control Charting under Drift,
in: \emph{Frontiers in Statistical Quality Control 10},
H.-J. Lenz, W. Schmid and P.-T. Wilrich (Eds.),
Physica Verlag, Heidelberg, Germany, 53-68.

C. Zou, Y. Liu and Z. Wang (2009),
Comparisons of control schemes for monitoring
the means of processes subject to drifts,
\emph{Metrika 70}, 141-163.
}
\author{Sven Knoth}
\seealso{
\code{xcusum.arl} and \code{xcusum.ad} for zero-state and
steady-state ARL computation of CUSUM control charts
for the classical step change model.
}
\examples{
## Gan (1992)
## Table 1
## original values are
#   deltas  arl
#   0.0001  475
#   0.0005  261
#   0.0010  187
#   0.0020  129
#   0.0050  76.3
#   0.0100  52.0
#   0.0200  35.2
#   0.0500  21.4
#   0.1000  15.0
#   0.5000  6.95
#   1.0000  5.16
#   3.0000  3.30
k <- .25
h <- 8
r <- 50
DxDcusum.arl <- Vectorize(xDcusum.arl, "delta")
deltas <- c(0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.5, 1, 3)
arl.like.Gan   <-
  round(DxDcusum.arl(k, h, deltas, r=r, with0=TRUE), digits=2)
arl.like.Knoth <-
  round(DxDcusum.arl(k, h, deltas, r=r, mode="Knoth", with0=TRUE), digits=2)
data.frame(deltas, arl.like.Gan, arl.like.Knoth)

## Zou et al. (2009)
## Table 1
## original values are
#  delta   arl1  arl2  arl3
#  0           ~ 1730
#  0.0005  345   412   470
#  0.001   231   275   317
#  0.005   86.6  98.6  112
#  0.01    56.9  61.8  69.3
#  0.05    22.6  21.6  22.7
#  0.1     15.4  14.7  14.2
#  0.5     6.60  5.54  5.17
#  1.0     4.63  3.80  3.45
#  2.0     3.17  2.67  2.32
#  3.0     2.79  2.04  1.96
#  4.0     2.10  1.98  1.74
\dontrun{
k1 <- 0.25
k2 <- 0.5
k3 <- 0.75
h1 <- 9.660
h2 <- 5.620
h3 <- 3.904
deltas <- c(0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1:4)
arl1 <- c(round(xcusum.arl(k1, h1, 0, r=r), digits=2),
          round(DxDcusum.arl(k1, h1, deltas, r=r), digits=2))
arl2 <- c(round(xcusum.arl(k2, h2, 0), digits=2),
          round(DxDcusum.arl(k2, h2, deltas, r=r), digits=2))
arl3 <- c(round(xcusum.arl(k3, h3, 0, r=r), digits=2),
          round(DxDcusum.arl(k3, h3, deltas, r=r), digits=2))
data.frame(delta=c(0, deltas), arl1, arl2, arl3)}
}
\keyword{ts}
