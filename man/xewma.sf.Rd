\name{xewma.sf}
\alias{xewma.sf}
\title{Compute the survival function of EWMA run length}
\description{Computation of the survival function of the Run Length (RL) for EWMA control charts monitoring normal mean.}
\usage{xewma.sf(l, c, mu, n, zr=0, hs=0, sided="one", limits="fix", q=1, r=40)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{mu}{true mean.}
\item{n}{calculate sf up to value \code{n}.}
\item{zr}{reflection border for the one-sided chart.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided EWMA control chart 
by choosing \code{"one"} and \code{"two"}, respectively.}
\item{limits}{distinguishes between different control limits behavior.}
\item{q}{change point position. For \eqn{q=1} and
\eqn{\mu=\mu_0} and \eqn{\mu=\mu_1}, the usual
zero-state situation for the in-control and out-of-control case, respectively,
are calculated. Note that mu0=0 is implicitely fixed.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-sided) or \code{r} (two-sided).}
}
\details{
The survival function P(L>n) and derived from it also the cdf P(L<=n) and the pmf P(L=n) illustrate
the distribution of the EWMA run length. For large n the geometric tail could be exploited. That is,
with reasonable large n the complete distribution is characterized.
The algorithm is based on Waldmann's survival function iteration procedure.
For varying limits and for change points after 1 the algorithm from Knoth (2004) is applied.
Note that for one-sided EWMA charts (\code{sided}=\code{"one"}), only
\code{"vacl"} and \code{"stat"} are deployed, while for two-sided ones
(\code{sided}=\code{"two"}) also \code{"fir"}, \code{"both"}
(combination of \code{"fir"} and \code{"vacl"}), and \code{"Steiner"} are implemented.
For details see Knoth (2004).
}
\value{Returns a vector which resembles the survival function up to a certain point.}
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
\code{xewma.arl} for zero-state ARL computation of EWMA control charts.
}
\examples{
## Gan (1993), two-sided EWMA with fixed control limits
## some values of his Table 1 -- any median RL should be 500

G.lambda <- c(.05, .1, .15, .2, .25)
G.h <- c(.441, .675, .863, 1.027, 1.177)/sqrt(G.lambda/(2-G.lambda))

for ( i in 1:length(G.lambda) ) {
  SF <- xewma.sf(G.lambda[i], G.h[i], 0, 1000)
  if (i==1) plot(1:length(SF), SF, type="l", xlab="n", ylab="P(L>n)")
  else lines(1:length(SF), SF, col=i)
  }
}
\keyword{ts}
