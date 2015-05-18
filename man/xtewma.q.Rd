\name{xtewma.q}
\alias{xtewma.q}
\alias{xtewma.q.crit}
\title{Compute RL quantiles of EWMA control charts}
\description{Computation of quantiles of the Run Length (RL)
for EWMA control charts monitoring normal mean.}
\usage{xtewma.q(l, c, df, mu, alpha, zr=0, hs=0, sided="two", limits="fix", mode="tan",
q=1, r=40)

xtewma.q.crit(l, L0, df, mu, alpha, zr=0, hs=0, sided="two", limits="fix", mode="tan",
r=40, c.error=1e-12, a.error=1e-9, OUTPUT=FALSE)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{df}{degrees of freedom -- parameter of the t distribution.}
\item{mu}{true mean.}
\item{alpha}{quantile level.}
\item{zr}{reflection border for the one-sided chart.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided EWMA control chart 
by choosing \code{"one"} and \code{"two"}, respectively.}
\item{limits}{distinguishes between different control limits behavior.}
\item{mode}{Controls the type of variables substitution that might improve the numerical performance. Currently,
\code{"identity"}, \code{"sin"}, \code{"sinh"}, and \code{"tan"} (default) are provided.}
\item{q}{change point position. For \eqn{q=1} and
\eqn{\mu=\mu_0} and \eqn{\mu=\mu_1}, the usual
zero-state ARLs for the in-control and out-of-control case, respectively,
are calculated. For \eqn{q>1} and \eqn{\mu!=0} conditional delays, that is,
\eqn{E_q(L-q+1|L\geq)}, will be determined.
Note that mu0=0 is implicitely fixed.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-sided) or \code{r}
(two-sided).}
\item{L0}{in-control quantile value.}
\item{c.error}{error bound for two succeeding values of the critical value
during applying the secant rule.}
\item{a.error}{error bound for the quantile level \code{alpha} during applying
the secant rule.}
\item{OUTPUT}{activate or deactivate additional output.}
}
\details{
Instead of the popular ARL (Average Run Length) quantiles of the EWMA
stopping time (Run Length) are determined. The algorithm is based on
Waldmann's survival function iteration procedure.
If \code{limits} is \code{"vacl"}, then the method presented in Knoth (2003) is utilized.
For details see Knoth (2004).
}
\value{Returns a single value which resembles the RL quantile of order \code{q}.}
\references{
F. F. Gan (1993),
An optimal design of EWMA control charts based on the median run length,
\emph{J. Stat. Comput. Simulation 45}, 169-184.

S. Knoth (2003),
EWMA schemes with non-homogeneous transition kernels,
\emph{Sequential Analysis 22}, 241-255.

S. Knoth (2004),
Fast initial response features for EWMA Control Charts,
\emph{Statistical Papers 46}, 47-64.

K.-H. Waldmann (1986),
Bounds for the distribution of the run length of geometric moving 
average charts, \emph{Appl. Statist. 35}, 151-158.
}
\author{Sven Knoth}
\seealso{
\code{xewma.q} for RL quantile computation of EWMA control charts in the normal case.
}
\examples{
## will follow
}
\keyword{ts}
