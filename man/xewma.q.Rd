\name{xewma.q}
\alias{xewma.q}
\alias{xewma.q.crit}
\title{Compute RL quantiles of EWMA control charts}
\description{Computation of quantiles of the Run Length (RL)
for EWMA control charts monitoring normal mean.}
\usage{xewma.q(l, c, mu, alpha, zr=0, hs=0, sided="two", limits="fix", q=1, r=40)

xewma.q.crit(l, L0, mu, alpha, zr=0, hs=0, sided="two", limits="fix", r=40,
c.error=1e-12, a.error=1e-9, OUTPUT=FALSE)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{mu}{true mean.}
\item{alpha}{quantile level.}
\item{zr}{reflection border for the one-sided chart.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided EWMA control chart 
by choosing \code{"one"} and \code{"two"}, respectively.}
\item{limits}{distinguishes between different control limits behavior.}
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
If \code{limits} is not \code{"fix"}, then the method presented
in Knoth (2003) is utilized.
Note that for one-sided EWMA charts (\code{sided}=\code{"one"}), only
\code{"vacl"} and \code{"stat"} are deployed, while for two-sided ones
(\code{sided}=\code{"two"}) also \code{"fir"}, \code{"both"}
(combination of \code{"fir"} and \code{"vacl"}), and \code{"Steiner"} are
implemented. For details see Knoth (2004).
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

S. Knoth (2015),
Run length quantiles of EWMA control charts monitoring normal mean or/and variance,
\emph{International Journal of Production Research 53}, 4629-4647.

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
XEWMA.Q <- Vectorize("xewma.q", c("l", "c"))
G.lambda <- c(.05, .1, .15, .2, .25)
G.h <- c(.441, .675, .863, 1.027, 1.177)
MEDIAN <- ceiling(XEWMA.Q(G.lambda, G.h/sqrt(G.lambda/(2-G.lambda)),
0, .5, sided="two"))
print(cbind(G.lambda, MEDIAN))

## increase accuracy of thresholds

# (i) calculate threshold for given in-control median value by
#     deplyoing secant rule
XEWMA.q.crit <- Vectorize("xewma.q.crit", "l")

# (ii) re-calculate the thresholds and remove the standardization step
L0 <- 500
G.h.new <- XEWMA.q.crit(G.lambda, L0, 0, .5, sided="two")
G.h.new <- round(G.h.new * sqrt(G.lambda/(2-G.lambda)), digits=5)

# (iii) compare Gan's original values and the new ones with 5 digits
print(cbind(G.lambda, G.h.new, G.h))

# (iv) calculate the new medians
MEDIAN <- ceiling(XEWMA.Q(G.lambda, G.h.new/sqrt(G.lambda/(2-G.lambda)),
0, .5, sided="two"))
print(cbind(G.lambda, MEDIAN))
}
\keyword{ts}
