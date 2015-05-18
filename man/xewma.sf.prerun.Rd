\name{xewma.sf.prerun}
\alias{xewma.sf.prerun}
\title{Compute the survival function of EWMA run length in case of estimated parameters}
\description{Computation of the survival function of the Run Length (RL) for EWMA
control charts monitoring normal mean
if the in-control mean, standard deviation, or both are estimated by a pre run.}
\usage{xewma.sf.prerun(l, c, mu, n, zr=0, hs=0, sided="one", limits="fix", q=1,
size=100, df=NULL, estimated="mu", qm.mu=30, qm.sigma=30,
truncate=1e-10, tail_approx=TRUE, bound=1e-10)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{mu}{true mean.}
\item{n}{calculate sf up to value \code{n}.}
\item{zr}{reflection border for the one-sided chart.}
\item{hs}{so-called headstart (give fast initial response).}
\item{sided}{distinguish between one- and two-sided EWMA control chart 
by choosing \code{"one"} and \code{"two"}, respectively.}
\item{limits}{distinguish between different control limits behavior.}
\item{q}{change point position. For \eqn{q=1} and
\eqn{\mu=\mu_0} and \eqn{\mu=\mu_1}, the usual
zero-state situation for the in-control and out-of-control case, respectively,
are calculated. Note that mu0=0 is implicitely fixed.}
\item{size}{pre run sample size.}
\item{df}{degrees of freedom of the pre run variance estimator. Typically it is simply \code{size} - 1.
If the pre run is collected in batches, then also other
values are needed.}
\item{estimated}{name the parameter to be estimated within the \code{"mu"},
\code{"sigma"}, \code{"both"}.}
\item{qm.mu}{number of quadrature nodes for convoluting the mean uncertainty.}
\item{qm.sigma}{number of quadrature nodes for convoluting the standard deviation uncertainty.}
\item{truncate}{size of truncated tail.}
\item{tail_approx}{Controls whether the geometric tail approximation is used (is faster) or not.}
\item{bound}{control when the geometric tail kicks in; the larger the quicker and less accurate; \code{bound} should be larger than 0 and less than 0.001.}
}
\details{
The survival function P(L>n) and derived from it also the cdf P(L<=n) and
the pmf P(L=n) illustrate the distribution of the EWMA run length...
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

L. A. Jones, C. W. Champ, S. E. Rigdon (2001),
The performance of exponentially weighted moving average charts
with estimated parameters,
\emph{Technometrics 43}, 156-167.

K.-H. Waldmann (1986),
Bounds for the distribution of the run length of geometric moving 
average charts, \emph{Appl. Statist. 35}, 151-158.
}
\author{Sven Knoth}
\seealso{
\code{xewma.sf} for the RL survival function of EWMA control charts
w/o pre run uncertainty.
}
\examples{
## Jones/Champ/Rigdon (2001)

c4m <- function(m, n) sqrt(2)*gamma( (m*(n-1)+1)/2 )/sqrt( m*(n-1) )/gamma( m*(n-1)/2 )

n <- 5 # sample size

# Figure 6, subfigure r=0.1
lambda <- 0.1
L <- 2.454

CDF0 <- 1 - xewma.sf(lambda, L, 0, 600, sided="two")
m <- 10 # pre run size
CDF1 <- 1 - xewma.sf.prerun(lambda, L/c4m(m,n), 0, 600, sided="two",
size=m, df=m*(n-1), estimated="both")
m <- 20
CDF2 <- 1 - xewma.sf.prerun(lambda, L/c4m(m,n), 0, 600, sided="two",
size=m, df=m*(n-1), estimated="both")
m <- 50
CDF3 <- 1 - xewma.sf.prerun(lambda, L/c4m(m,n), 0, 600, sided="two",
size=m, df=m*(n-1), estimated="both")

plot(CDF0, type="l", xlab="t", ylab=expression(P(T<=t)), xlim=c(0,500), ylim=c(0,1))
abline(v=0, h=c(0,1), col="grey", lwd=.7)
points((1:5)*100, CDF0[(1:5)*100], pch=18)
lines(CDF1, col="blue")
points((1:5)*100, CDF1[(1:5)*100], pch=2, col="blue")
lines(CDF2, col="red")
points((1:5)*100, CDF2[(1:5)*100], pch=16, col="red")
lines(CDF3, col="green")
points((1:5)*100, CDF3[(1:5)*100], pch=5, col="green")

legend("bottomright", c("Known", "m=10, n=5", "m=20, n=5", "m=50, n=5"),
       col=c("black", "blue", "red", "green"), pch=c(18, 2, 16, 5), lty=1)
}
\keyword{ts}
