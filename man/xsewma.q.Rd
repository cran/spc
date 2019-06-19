\name{xsewma.q}
\alias{xsewma.q}
\alias{xsewma.q.crit}
\title{Compute critical values of simultaneous EWMA control charts
(mean and variance charts) for given RL quantile}
\description{Computation of the critical values (similar to alarm limits)
for different types of simultaneous EWMA control charts 
(based on the sample mean and the sample variance \eqn{S^2})
monitoring normal mean and variance.}
\usage{xsewma.q(lx, cx, ls, csu, df, alpha, mu, sigma, hsx=0,
Nx=40, csl=0, hss=1, Ns=40, sided="upper", qm=30)

xsewma.q.crit(lx, ls, L0, alpha, df, mu0=0, sigma0=1, csu=NULL,
hsx=0, hss=1, sided="upper", mode="fixed", Nx=20, Ns=40, qm=30,
c.error=1e-12, a.error=1e-9)}
\arguments{
\item{lx}{smoothing parameter lambda of the two-sided mean EWMA chart.}
\item{cx}{control limit of the two-sided mean EWMA control chart.}
\item{ls}{smoothing parameter lambda of the variance EWMA chart.}
\item{csu}{for two-sided (\code{sided}=\code{"two"}) and fixed upper
control limit (\code{mode}=\code{"fixed"}, only for \code{xsewma.q.crit})
a value larger than \code{sigma0}
has to been given, for all other cases \code{cu} is ignored.
It is the upper control limit of the variance EWMA control chart.}
\item{L0}{in-control RL quantile at level \code{alpha}.}
\item{df}{actual degrees of freedom, corresponds to subgroup size
(for known mean it is equal to the subgroup size,
for unknown mean it is equal to subgroup size minus one.}
\item{alpha}{quantile level.}
\item{mu}{true mean.}
\item{sigma}{true standard deviation.}
\item{mu0}{in-control mean.}
\item{sigma0}{in-control standard deviation.}
\item{hsx}{so-called headstart (enables fast initial response) of the mean chart -- 
do not confuse with the true FIR feature considered in xewma.arl; will be updated.}
\item{Nx}{dimension of the approximating matrix of the mean chart.}
\item{csl}{lower control limit of the variance EWMA control chart; default value is 0;
not considered if \code{sided} is \code{"upper"}.}
\item{hss}{headstart (enables fast initial response) of the variance chart.}
\item{Ns}{dimension of the approximating matrix of the variance chart.}
\item{sided}{distinguishes between one- and two-sided two-sided
EWMA-\eqn{S^2}{S^2} control charts by choosing \code{"upper"} (upper chart without
reflection at \code{cl} -- the actual value of of \code{cl} is not used).}
\item{mode}{only deployed for \code{sided}=\code{"two"} -- with \code{"fixed"}
an upper control limit (see \code{cu}) is set and only the lower is
determined to obtain the in-control ARL \code{L0}, while with \code{"unbiased"}
a certain unbiasedness of the ARL function is guaranteed (here, both the
lower and the upper control limit are calculated).}
\item{qm}{number of quadrature nodes used for the collocation integrals.}
\item{c.error}{error bound for two succeeding values of the critical value
during applying the secant rule.}
\item{a.error}{error bound for the quantile level \code{alpha} during
applying the secant rule.}
}
\details{
Instead of the popular ARL (Average Run Length) quantiles of the EWMA
stopping time (Run Length) are determined. The algorithm is based on
Waldmann's survival function iteration procedure and on Knoth (2007).
\code{xsewma.q.crit} determines the critical values (similar to alarm limits)
for given in-control RL quantile \code{L0} at level \code{alpha} by applying secant
rule and using \code{xsewma.sf()}.
In case of \code{sided}=\code{"two"} and \code{mode}=\code{"unbiased"}
a two-dimensional secant rule is applied that also ensures that the
maximum of the RL cdf for given standard deviation is attained at \code{sigma0}.
}
\value{Returns a single value which resembles the RL quantile of order \code{alpha} and
the critical value of the two-sided mean EWMA chart and
the lower and upper controls limit \code{csl} and \code{csu} of the
variance EWMA chart, respectively.}
\references{

S. Knoth (2007),
Accurate ARL calculation for EWMA control charts monitoring simultaneously
normal mean and variance,
\emph{Sequential Analysis 26}, 251-264.
}
\author{Sven Knoth}
\seealso{\code{xsewma.arl} for calculation of ARL of simultaneous EWMA charts and
\code{xsewma.sf} for the RL survival function.}
\examples{
## will follow
}
\keyword{ts}
