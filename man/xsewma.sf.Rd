\name{xsewma.sf}
\alias{xsewma.sf}
\title{Compute the survival function of simultaneous EWMA control
charts (mean and variance charts)}
\description{Computation of the survival function of the Run Length (RL)
for EWMA control charts monitoring simultaneously normal mean and variance.}
\usage{xsewma.sf(n, lx, cx, ls, csu, df, mu, sigma, hsx=0, Nx=40,
csl=0, hss=1, Ns=40, sided="upper", qm=30) }
\arguments{
\item{n}{calculate sf up to value \code{n}.}
\item{lx}{smoothing parameter lambda of the two-sided mean EWMA chart.}
\item{cx}{control limit of the two-sided mean EWMA control chart.}
\item{ls}{smoothing parameter lambda of the variance EWMA chart.}
\item{csu}{upper control limit of the variance EWMA control chart.}
\item{df}{actual degrees of freedom, corresponds to subgroup size (for known mean it is
equal to the subgroup size, for unknown mean it is equal to subgroup size minus one.}
\item{mu}{true mean.}
\item{sigma}{true standard deviation.}
\item{hsx}{so-called headstart (enables fast initial response) of the mean chart --
do not confuse with the true FIR feature considered in xewma.arl; will be updated.}
\item{Nx}{dimension of the approximating matrix of the mean chart.}
\item{csl}{lower control limit of the variance EWMA control chart; default value is 0;
not considered if \code{sided} is \code{"upper"}.}
\item{hss}{headstart (enables fast initial response) of the variance chart.}
\item{Ns}{dimension of the approximating matrix of the variance chart.}
\item{sided}{distinguishes between one- and two-sided two-sided
EWMA-\eqn{S^2}{S^2} control charts by choosing \code{"upper"} (upper chart
without reflection at \code{cl} -- the actual value of
\code{cl} is not used), \code{"Rupper"} (upper chart with reflection at \code{cl}),
\code{"Rlower"} (lower chart with reflection at \code{cu}), and \code{"two"}
(two-sided chart), respectively.} 
\item{qm}{number of quadrature nodes used for the collocation integrals.}
}
\details{
The survival function P(L>n) and derived from it also the cdf P(L<=n) and
the pmf P(L=n) illustrate
the distribution of the EWMA run length. For large n the geometric tail
could be exploited. That is,
with reasonable large n the complete distribution is characterized.
The algorithm is based on Waldmann's survival function iteration procedure and
on results in Knoth (2007).
}
\value{Returns a vector which resembles the survival function up to a certain point.}
\references{
S. Knoth (2007),
Accurate ARL calculation for EWMA control charts monitoring simultaneously
normal mean and variance,
\emph{Sequential Analysis 26}, 251-264.

K.-H. Waldmann (1986),
Bounds for the distribution of the run length of geometric moving 
average charts, \emph{Appl. Statist. 35}, 151-158.
}
\author{Sven Knoth}
\seealso{
\code{xsewma.arl} for zero-state ARL computation of simultaneous EWMA
control charts.
}
\examples{
## will follow
}
\keyword{ts}
