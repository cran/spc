\name{xDgrsr.arl}
\alias{xDgrsr.arl}
\title{Compute ARLs of Shiryaev-Roberts schemes under drift}
\description{Computation of the (zero-state and other)
Average Run Length (ARL) under drift
for Shiryaev-Roberts schemes monitoring normal mean.}
\usage{xDgrsr.arl(k, g, delta, zr = 0, hs = NULL, sided = "one", m = NULL,
mode = "Gan", q = 1, r = 30, with0 = FALSE)}
\arguments{
\item{k}{reference value of the Shiryaev-Roberts scheme.}
\item{g}{control limit (alarm threshold) of Shiryaev-Roberts scheme.}
\item{delta}{true drift parameter.}
\item{zr}{reflection border for the one-sided chart.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided 
Shiryaev-Roberts schemes
by choosing \code{"one"} and \code{"two"}, respectively. Currentlly,
the two-sided scheme is not implemented.}
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
\item{mode}{decide whether Gan's or Knoth's approach is used. Use
\code{"Gan"} and \code{"Knoth"}, respectively.
\code{"Knoth"} is not implemented yet.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-sided) or \code{r} (two-sided).}
\item{with0}{defines whether the first observation used for the RL calculation
follows already 1*delta or still 0*delta.
With \code{q} additional flexibility is given.}
}
\details{
Based on Gan (1991) or Knoth (2003), the ARL is calculated for Shiryaev-Roberts schemes under drift.
In case of Gan's framework, the usual
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
\code{xewma.arl} and \code{xewma.ad} for zero-state and
steady-state ARL computation of EWMA control charts
for the classical step change model.
}
\examples{
\dontrun{
## Monte Carlo example with 10^8 replicates
#   delta      arl    s.e.
#   0.0001 381.8240   0.0304
#   0.0005 238.4630   0.0148
#   0.001  177.4061   0.0097
#   0.002  125.9055   0.0061
#   0.005   75.7574   0.0031
#   0.01    50.2203   0.0018
#   0.02    32.9458   0.0011
#   0.05    18.9213   0.0005
#   0.1     12.6054   0.0003
#   0.5      5.2157   0.0001
#   1        3.6537   0.0001
#   3        2.0289   0.0000
k <- .5
L0 <- 500
zr <- -7
r <- 50
g <- xgrsr.crit(k, L0, zr=zr, r=r)
DxDgrsr.arl <- Vectorize(xDgrsr.arl, "delta")
deltas <- c(0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.5, 1, 3)
arls <- round(DxDgrsr.arl(k, g, deltas, zr=zr, r=r), digits=4)
data.frame(deltas, arls)
}
}
\keyword{ts}
