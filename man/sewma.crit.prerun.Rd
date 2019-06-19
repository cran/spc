\name{sewma.crit.prerun}
\alias{sewma.crit.prerun}
\title{Compute critical values of of EWMA (variance charts) control charts under pre-run uncertainty}
\description{Computation of quantiles of the Run Length (RL) for EWMA control
charts monitoring normal variance.}
\usage{sewma.crit.prerun(l,L0,df1,df2,sigma0=1,cl=NULL,cu=NULL,hs=1,sided="upper",
mode="fixed",r=40,qm=30,qm.sigma=30,truncate=1e-10,
tail_approx=TRUE,c.error=1e-10,a.error=1e-9)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{L0}{in-control quantile value.}
\item{df1}{actual degrees of freedom, corresponds to subgroup size
(for known mean it is equal to the subgroup size,
for unknown mean it is equal to subgroup size minus one.}
\item{df2}{degrees of freedom of the pre-run variance estimator.}
\item{sigma,sigma0}{true and in-control standard deviation, respectively.}
\item{cl}{deployed for \code{sided}=\code{"Rupper"}, that is, upper variance control chart with lower
reflecting barrier \code{cl}.}
\item{cu}{for two-sided (\code{sided}=\code{"two"}) and fixed upper control limit
(\code{mode}=\code{"fixed"}) a value larger than \code{sigma0}
has to been given, for all other cases \code{cu} is ignored.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided two-sided EWMA-\eqn{S^2}{S^2} control charts
by choosing \code{"upper"} (upper chart without reflection at \code{cl} -- the actual value of \code{cl}
is not used), \code{"Rupper"} (upper chart with reflection at \code{cl}), \code{"Rlower"} (lower chart
with reflection at \code{cu}),and \code{"two"} (two-sided chart), respectively.}
\item{mode}{only deployed for \code{sided}=\code{"two"} -- with \code{"fixed"} an upper control limit
(see \code{cu}) is set and only the lower is calculated to obtain the in-control ARL \code{L0}, while
with \code{"unbiased"} a certain unbiasedness of the ARL function is guaranteed (here, both the
lower and the upper control limit are calculated).}
\item{r}{dimension of the resulting linear equation system (highest order
of the collocation polynomials).}
\item{qm}{number of quadrature nodes for calculating the collocation definite integrals.}
\item{qm.sigma}{number of quadrature nodes for convoluting the standard deviation uncertainty.}
\item{truncate}{size of truncated tail.}
\item{tail_approx}{controls whether the geometric tail approximation is used (is faster) or not.}
\item{c.error}{error bound for two succeeding values of the critical value
during applying the secant rule.}
\item{a.error}{error bound for the quantile level \code{alpha} during applying the secant rule.}
}
\details{
\code{sewma.crit.prerun} determines the critical values (similar to alarm limits)
for given in-control ARL \code{L0}
by applying secant rule and using \code{sewma.arl.prerun()}.
In case of \code{sided}=\code{"two"} and \code{mode}=\code{"unbiased"}
a two-dimensional secant rule is applied that also ensures that the
maximum of the ARL function for given standard deviation is attained
at \code{sigma0}. See Knoth (2010) for some details of the algorithm involved.
}
\value{Returns the lower and upper control limit \code{cl} and \code{cu}.}
\references{
H.-J. Mittag and D. Stemann and B. Tewes (1998),
EWMA-Karten zur \"Uberwachung der Streuung von Qualit\"atsmerkmalen,
\emph{Allgemeines Statistisches Archiv 82}, 327-338,
S. Knoth (2005),
Accurate ARL computation for EWMA-\eqn{S^2}{S^2} control charts,
\emph{Statistics and Computing 15}, 341-352.

S. Knoth (2010),
Control Charting Normal Variance -- Reflections, Curiosities, and Recommendations,
in \emph{Frontiers in Statistical Quality Control 9},
H.-J. Lenz and P.-T. Wilrich (Eds.),
Physica Verlag, Heidelberg, Germany, 3-18.
}
\author{Sven Knoth}
\seealso{\code{sewma.arl.prerun} for calculation of ARL of variance charts under
pre-run uncertainty and \code{sewma.crit} for
the algorithm w/o pre-run uncertainty.}
\examples{
## will follow
}
\keyword{ts}
