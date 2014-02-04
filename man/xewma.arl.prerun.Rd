\name{xewma.arl.prerun}
\alias{xewma.arl.prerun}
\alias{xewma.crit.prerun}
\title{Compute ARLs of EWMA control charts in case of estimated parameters}
\description{Computation of the (zero-state) Average Run Length (ARL)
for different types of EWMA control charts monitoring normal mean
if the in-control mean, standard deviation, or both are estimated by a pre run.}
\usage{xewma.arl.prerun(l, c, mu, zr=0, hs=0, sided="two", limits="fix", q=1,
size=100, df=NULL, estimated="mu", qm.mu=30, qm.sigma=30, truncate=1e-10)

xewma.crit.prerun(l, L0, mu, zr=0, hs=0, sided="two", limits="fix", size=100,
df=NULL, estimated="mu", qm.mu=30, qm.sigma=30, truncate=1e-10,
c.error=1e-12, L.error=1e-9, OUTPUT=FALSE)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{mu}{true mean shift.}
\item{zr}{reflection border for the one-sided chart.}
\item{hs}{so-called headstart (give fast initial response).}
\item{sided}{distinguish between one- and two-sided EWMA control chart 
by choosing \code{"one"} and \code{"two"}, respectively.}
\item{limits}{distinguish between different control limits behavior.}
\item{q}{change point position. For \eqn{q=1} and
\eqn{\mu=\mu_0} and \eqn{\mu=\mu_1}, the usual
zero-state ARLs for the in-control and out-of-control case, respectively,
are calculated. For \eqn{q>1} and \eqn{\mu!=0} conditional delays, that is,
\eqn{E_q(L-q+1|L\ge q)}, will be determined.
Note that mu0=0 is implicitely fixed.}
\item{size}{pre run sample size.}
\item{df}{Degrees of freedom of the pre run variance estimator. Typically it is simply \code{size} - 1.
If the pre run is collected in batches, then also other values are needed.}
\item{estimated}{name the parameter to be estimated within
the \code{"mu"}, \code{"sigma"}, \code{"both"}.}
\item{qm.mu}{number of quadrature nodes for convoluting the mean uncertainty.}
\item{qm.sigma}{number of quadrature nodes for convoluting the standard deviation uncertainty.}
\item{truncate}{size of truncated tail.}
\item{L0}{in-control ARL.}
\item{c.error}{error bound for two succeeding values of the critical value
during applying the secant rule.}
\item{L.error}{error bound for the ARL level \code{L0} during applying the secant rule.}
\item{OUTPUT}{activate or deactivate additional output.}
}
\details{
Essentially, the ARL function \code{xewma.arl} is convoluted with the
distribution of the sample mean, standard deviation or both.
For details see Jones/Champ/Rigdon (2001) and Knoth (2014?).
}
\value{Returns a single value which resembles the ARL.}
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
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} for the usual zero-state ARL computation.
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

xewma.ARL <- Vectorize("xewma.arl", "mu")
xewma.ARL.prerun <- Vectorize("xewma.arl.prerun", "mu")

mu <- c(0, .25, .5, 1, 1.5, 2)
ARL <- round(xewma.ARL(lambda, L, mu, sided="two"), digits=2)
p.ARL <- round(xewma.ARL.prerun(lambda, L/C4m, mu, sided="two",
size=m, df=m*(n-1), estimated="both", qm.mu=70), digits=2)

# Monte-Carlo with 10^8 repetitions: 200.325 (0.020) and 144.458 (0.022)
cbind(mu, ARL, p.ARL)

\dontrun{
# Figure 5, subfigure r = 0.2
mu_ <- (0:85)/40
ARL_ <- round(xewma.ARL(lambda, L, mu_, sided="two"), digits=2)
p.ARL_ <- round(xewma.ARL.prerun(lambda, L/C4m, mu_, sided="two",
size=m, df=m*(n-1), estimated="both"), digits=2)

plot(mu_, ARL_, type="l", xlab=expression(delta), ylab="ARL", xlim=c(0,2))
abline(v=0, h=0, col="grey", lwd=.7)
points(mu, ARL, pch=5)
lines(mu_, p.ARL_, col="blue")
points(mu, p.ARL, pch=18, col ="blue")
legend("topright", c("Known", "Estimated"), col=c("black", "blue"),
lty=1, pch=c(5, 18))
}
}
\keyword{ts}
