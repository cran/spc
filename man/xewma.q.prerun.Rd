\name{xewma.q.prerun}
\alias{xewma.q.prerun}
\alias{xewma.q.crit.prerun}
\title{Compute RL quantiles of EWMA control charts in case of estimated parameters}
\description{Computation of quantiles of the Run Length (RL)
for EWMA control charts monitoring normal mean
if the in-control mean, standard deviation, or both are estimated by a pre run.}
\usage{xewma.q.prerun(l, c, mu, p, zr=0, hs=0, sided="two", limits="fix", q=1, size=100,
df=NULL, estimated="mu", qm.mu=30, qm.sigma=30, truncate=1e-10, bound=1e-10)

xewma.q.crit.prerun(l, L0, mu, p, zr=0, hs=0, sided="two", limits="fix", size=100,
df=NULL, estimated="mu", qm.mu=30, qm.sigma=30, truncate=1e-10, bound=1e-10,
c.error=1e-10, p.error=1e-9, OUTPUT=FALSE)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{mu}{true mean shift.}
\item{p}{quantile level.}
\item{zr}{reflection border for the one-sided chart.}
\item{hs}{so-called headstart (give fast initial response).}
\item{sided}{distinguish between one- and two-sided EWMA control chart 
by choosing \code{"one"} and \code{"two"}, respectively.}
\item{limits}{distinguish between different control limits behavior.}
\item{q}{change point position. For \eqn{q=1} and
\eqn{\mu=\mu_0} and \eqn{\mu=\mu_1}, the usual
zero-state ARLs for the in-control and out-of-control case, respectively,
are calculated. For \eqn{q>1} and \eqn{\mu!=0} conditional delays, that is,
\eqn{E_q(L-q+1|L\geq)}, will be determined.
Note that mu0=0 is implicitely fixed.}
\item{size}{pre run sample size.}
\item{df}{Degrees of freedom of the pre run variance estimator. Typically it is simply \code{size} - 1.
If the pre run is collected in batches, then also other values are needed.}
\item{estimated}{name the parameter to be estimated within the \code{"mu"}, \code{"sigma"},
\code{"both"}.}
\item{qm.mu}{number of quadrature nodes for convoluting the mean uncertainty.}
\item{qm.sigma}{number of quadrature nodes for convoluting the standard deviation uncertainty.}
\item{truncate}{size of truncated tail.}
\item{bound}{control when the geometric tail kicks in; the larger the quicker and less accurate; \code{bound} should be larger than 0 and less than 0.001.}
\item{L0}{in-control quantile value.}
\item{c.error}{error bound for two succeeding values of the critical value during
applying the secant rule.}
\item{p.error}{error bound for the quantile level \code{p} during applying the secant rule.}
\item{OUTPUT}{activate or deactivate additional output.}
}
\details{
Essentially, the ARL function \code{xewma.q} is convoluted with the
distribution of the sample mean, standard deviation or both.
For details see Jones/Champ/Rigdon (2001) and Knoth (2014?).
}
\value{Returns a single value which resembles the RL quantile of order \code{q}.}
\references{
L. A. Jones, C. W. Champ, S. E. Rigdon (2001),
The performance of exponentially weighted moving average charts
with estimated parameters,
\emph{Technometrics 43}, 156-167.

S. Knoth (2003),
EWMA schemes with non-homogeneous transition kernels,
\emph{Sequential Analysis 22}, 241-255.

S. Knoth (2004),
Fast initial response features for EWMA Control Charts,
\emph{Statistical Papers 46}, 47-64.

S. Knoth (2014?),
tbd,
\emph{tbd}, tbd-tbd.

K.-H. Waldmann (1986),
Bounds for the distribution of the run length of geometric moving 
average charts, \emph{Appl. Statist. 35}, 151-158.
}
\author{Sven Knoth}
\seealso{
\code{xewma.q} for the usual RL quantiles computation of EWMA control charts.
}
\examples{
## Jones/Champ/Rigdon (2001)

c4m <- function(m, n) sqrt(2)*gamma( (m*(n-1)+1)/2 )/sqrt( m*(n-1) )/gamma( m*(n-1)/2 )

n <- 5 # sample size
m <- 20 # pre run with 20 samples of size n = 5
C4m <- c4m(m, n) # needed for bias correction

# Table 1, 3rd column
lambda <- 0.2
L <- 2.636

xewma.Q <- Vectorize("xewma.q", "mu")
xewma.Q.prerun <- Vectorize("xewma.q.prerun", "mu")

mu <- c(0, .25, .5, 1, 1.5, 2)
Q1  <- ceiling(xewma.Q(lambda, L, mu, 0.1, sided="two"))
Q2  <- ceiling(xewma.Q(lambda, L, mu, 0.5, sided="two"))
Q3  <- ceiling(xewma.Q(lambda, L, mu, 0.9, sided="two"))

cbind(mu, Q1, Q2, Q3)

\dontrun{
p.Q1 <- xewma.Q.prerun(lambda, L/C4m, mu, 0.1, sided="two", 
size=m, df=m*(n-1), estimated="both")
p.Q2 <- xewma.Q.prerun(lambda, L/C4m, mu, 0.5, sided="two",
size=m, df=m*(n-1), estimated="both")
p.Q3 <- xewma.Q.prerun(lambda, L/C4m, mu, 0.9, sided="two",
size=m, df=m*(n-1), estimated="both")

cbind(mu, p.Q1, p.Q2, p.Q3)
}

## original values are
#    mu Q1  Q2  Q3 p.Q1 p.Q2 p.Q3
#  0.00 25 140 456   13   73  345
#  0.25 12  56 174    9   46  253
#  0.50  7  20  56    6   20  101
#  1.00  4   7  15    3    7   18
#  1.50  3   4   7    2    4    8
#  2.00  2   3   5    2    3    5
}
\keyword{ts}
