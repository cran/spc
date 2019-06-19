\name{sewma.sf}
\alias{sewma.sf}
\title{Compute the survival function of EWMA run length}
\description{Computation of the survival function of the Run Length (RL)
for EWMA control charts monitoring normal variance.}
\usage{sewma.sf(n, l, cl, cu, sigma, df, hs=1, sided="upper", r=40, qm=30)}
\arguments{
\item{n}{calculate sf up to value \code{n}.}
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{cl}{lower control limit of the EWMA control chart.}
\item{cu}{upper control limit of the EWMA control chart.}
\item{sigma}{true standard deviation.}
\item{df}{actual degrees of freedom, corresponds to subgroup size
(for known mean it is equal to the subgroup size,
for unknown mean it is equal to subgroup size minus one.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided two-sided
EWMA-\eqn{S^2}{S^2} control charts 
by choosing \code{"upper"} (upper chart without reflection at \code{cl} --
the actual value of \code{cl} is not used),
\code{"Rupper"} (upper chart with reflection at \code{cl}), \code{"Rlower"}
(lower chart with reflection at \code{cu}),
and \code{"two"} (two-sided chart), respectively.} 
\item{r}{dimension of the resulting linear equation system
(highest order of the collocation polynomials).}
\item{qm}{number of quadrature nodes for calculating the
collocation definite integrals.}
}
\details{
The survival function P(L>n) and derived from it also the cdf P(L<=n)
and the pmf P(L=n) illustrate
the distribution of the EWMA run length. For large n the
geometric tail could be exploited. That is,
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
\code{sewma.arl} for zero-state ARL computation of variance EWMA control charts.
}
\examples{
## will follow
}
\keyword{ts}
