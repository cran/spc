\name{xcusum.q}
\alias{xcusum.q}
\title{Compute RL quantiles of CUSUM control charts}
\description{Computation of quantiles of the Run Length (RL)for CUSUM control charts monitoring normal mean.}
\usage{xcusum.q(k, h, mu, p, hs=0, sided="one", r=40)}
\arguments{
\item{k}{reference value of the CUSUM control chart.}
\item{h}{decision interval (alarm limit, threshold) of the CUSUM control chart.}
\item{mu}{true mean.}
\item{p}{quantile level.}
\item{hs}{so-called headstart (give fast initial response).}
\item{sided}{distinguish between one- and two-sided CUSUM control chart by choosing \code{"one"} and \code{"two"}, respectively.}
\item{r}{number of quadrature nodes, dimension of the resulting linear equation system is equal to \code{r+1}.}
}
\details{
Instead of the popular ARL (Average Run Length) quantiles of the CUSUM
stopping time (Run Length) are determined. The algorithm is based on
Waldmann's survival function iteration procedure.
}
\value{Returns a single value which resembles the RL quantile of order \code{q}.}
\references{
K.-H. Waldmann (1986),
Bounds for the distribution of the run length of one-sided and two-sided CUSUM quality control schemes,
\emph{Technometrics 28}, 61-67.
}
\author{Sven Knoth}
\seealso{
\code{xcusum.arl} for zero-state ARL computation of CUSUM control charts.
}
\examples{
## Waldmann (1986), one-sided CUSUM, Table 2
## original values are 345, 82, 9

XCUSUM.Q <- Vectorize("xcusum.q", "p")
k <- .5
h <- 3
mu <- 0 # corresponds to Waldmann's -0.5
p.list <- c(.95, .5, .05)
rl.quantiles <- ceiling(XCUSUM.Q(k, h, mu, p.list))
cbind(p.list, rl.quantiles)
}
\keyword{ts}
