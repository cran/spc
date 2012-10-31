\name{xcusum.sf}
\alias{xcusum.sf}
\title{Compute the survival function of CUSUM run length}
\description{Computation of the survival function of the Run Length (RL) for CUSUM control charts monitoring normal mean.}
\usage{xcusum.sf(k, h, mu, n, hs=0, sided="one", r=40)}
\arguments{
\item{k}{reference value of the CUSUM control chart.}
\item{h}{decision interval (alarm limit, threshold) of the CUSUM control chart.}
\item{mu}{true mean.}
\item{n}{calculate sf up to value \code{n}.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided CUSUM control chart by choosing \code{"one"} and \code{"two"}, respectively.}
\item{r}{number of quadrature nodes, dimension of the resulting linear equation system is equal to \code{r+1}.}
}
\details{
The survival function P(L>n) and derived from it also the cdf P(L<=n) and the pmf P(L=n) illustrate
the distribution of the CUSUM run length. For large n the geometric tail could be exploited. That is,
with reasonable large n the complete distribution is characterized.
The algorithm is based on Waldmann's survival function iteration procedure.
}
\value{Returns a vector which resembles the survival function up to a certain point.}
\references{
K.-H. Waldmann (1986),
Bounds for the distribution of the run length of one-sided and two-sided CUSUM quality control schemes,
\emph{Technometrics 28}, 61-67.
}
\author{Sven Knoth}
\seealso{
\code{xcusum.q} for computation of CUSUM run length quantiles.
}
\examples{
## Waldmann (1986), one-sided CUSUM, Table 2

k <- .5
h <- 3
mu <- 0 # corresponds to Waldmann's -0.5
SF <- xcusum.sf(k, h, 0, 1000)
plot(1:length(SF), SF, type="l", xlab="n", ylab="P(L>n)", ylim=c(0,1))
#
}
\keyword{ts}
