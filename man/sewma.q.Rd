\name{sewma.q}
\alias{sewma.q}
\alias{sewma.q.crit}
\title{Compute RL quantiles of EWMA (variance charts) control charts}
\description{Computation of quantiles of the Run Length (RL)
for EWMA control charts monitoring normal variance.}
\usage{sewma.q(l, cl, cu, sigma, df, alpha, hs=1, sided="upper", r=40, qm=30)

sewma.q.crit(l,L0,alpha,df,sigma0=1,cl=NULL,cu=NULL,hs=1,sided="upper",
mode="fixed",ur=4,r=40,qm=30,c.error=1e-12,a.error=1e-9)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{cl}{deployed for \code{sided}=\code{"Rupper"}, that is, upper variance control
chart with lower reflecting barrier \code{cl}.}
\item{cu}{for two-sided (\code{sided}=\code{"two"}) and fixed upper control limit
(\code{mode}=\code{"fixed"}) a value larger than \code{sigma0}
has to been given, for all other cases \code{cu} is ignored.}
\item{sigma,sigma0}{true and in-control standard deviation, respectively.}
\item{df}{actual degrees of freedom, corresponds to subgroup size (for known mean it
is equal to the subgroup size,
for unknown mean it is equal to subgroup size minus one.}
\item{alpha}{quantile level.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided two-sided EWMA-\eqn{S^2}{S^2}
control charts by choosing \code{"upper"} (upper chart
without reflection at \code{cl} -- the actual value of \code{cl} is not used), \code{"Rupper"}
(upper chart with reflection at \code{cl}),
\code{"Rlower"} (lower chart with reflection at \code{cu}),and \code{"two"}
(two-sided chart), respectively.}
\item{mode}{only deployed for \code{sided}=\code{"two"} -- with \code{"fixed"} an upper
control limit (see \code{cu}) is set and only the lower is
calculated to obtain the in-control ARL \code{L0}, while with \code{"unbiased"} a
certain unbiasedness of the ARL function is guaranteed (here, both the
lower and the upper control limit are calculated).}
\item{ur}{truncation of lower chart for \code{classic} mode.}
\item{r}{dimension of the resulting linear equation system (highest order of the
collocation polynomials).}
\item{qm}{number of quadrature nodes for calculating the collocation definite integrals.}
\item{L0}{in-control quantile value.}
\item{c.error}{error bound for two succeeding values of the critical value during
applying the secant rule.}
\item{a.error}{error bound for the quantile level \code{alpha}
during applying the secant rule.}
}
\details{
Instead of the popular ARL (Average Run Length) quantiles of the EWMA
stopping time (Run Length) are determined. The algorithm is based on
Waldmann's survival function iteration procedure.
Thereby the ideas presented in Knoth (2007) are used.
\code{sewma.q.crit} determines the critical values (similar to alarm limits)
for given in-control RL quantile \code{L0} at level \code{alpha} by applying
secant rule and using \code{sewma.sf()}.
In case of \code{sided}=\code{"two"} and \code{mode}=\code{"unbiased"} a two-dimensional
secant rule is applied that also ensures that the
minimum of the cdf for given standard deviation is attained at \code{sigma0}.
}
\value{Returns a single value which resembles the RL quantile of order \code{alpha} and
the lower and upper control limit \code{cl} and \code{cu}, respectively.}
\references{
H.-J. Mittag and D. Stemann and B. Tewes (1998),
EWMA-Karten zur \"Uberwachung der Streuung von Qualit\"atsmerkmalen,
\emph{Allgemeines Statistisches Archiv 82}, 327-338,

C. A. Acosta-Mej\'ia and J. J. Pignatiello Jr. and B. V. Rao (1999),
A comparison of control charting procedures for monitoring process dispersion,
\emph{IIE Transactions 31}, 569-579.

S. Knoth (2005),
Accurate ARL computation for EWMA-\eqn{S^2}{S^2} control charts,
\emph{Statistics and Computing 15}, 341-352.

S. Knoth (2007),
Accurate ARL calculation for EWMA control charts monitoring simultaneously normal mean and variance,
\emph{Sequential Analysis 26}, 251-264.

S. Knoth (2010),
Control Charting Normal Variance -- Reflections, Curiosities, and Recommendations,
in \emph{Frontiers in Statistical Quality Control 9},
H.-J. Lenz and P.-T. Wilrich (Eds.),
Physica Verlag, Heidelberg, Germany, 3-18.

K.-H. Waldmann (1986),
Bounds for the distribution of the run length of geometric moving 
average charts, \emph{Appl. Statist. 35}, 151-158.
}
\author{Sven Knoth}
\seealso{\code{sewma.arl} for calculation of ARL of variance charts and
\code{sewma.sf} for the RL survival function.}
\examples{
## will follow
}
\keyword{ts}
