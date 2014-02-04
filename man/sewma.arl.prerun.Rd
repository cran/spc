\name{sewma.arl.prerun}
\alias{sewma.arl.prerun}
\title{Compute ARLs of EWMA control charts (variance charts) in case of estimated parameters}
\description{Computation of the (zero-state) Average Run Length (ARL)
for EWMA control charts  (based on the sample variance \eqn{S^2})
monitoring normal variance with estimated parameters.}
\usage{sewma.arl.prerun(l, cl, cu, sigma, df1, df2, hs=1, sided="upper",
r=40, qm=30, qm.sigma=30, truncate=1e-10)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{cl}{lower control limit of the EWMA control chart.}
\item{cu}{upper control limit of the EWMA control chart.}
\item{sigma}{true standard deviation.}
\item{df1}{actual degrees of freedom, corresponds to subgroup size (for known mean it is equal to the subgroup size,
for unknown mean it is equal to subgroup size minus one.}
\item{df2}{degrees of freedom of the pre-run variance estimator.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided two-sided EWMA-\eqn{S^2}{S^2} control charts 
by choosing \code{"upper"} (upper chart without reflection at \code{cl} -- the actual value of
\code{cl} is not used), \code{"Rupper"} (upper chart with reflection at \code{cl}),\code{"Rlower"}
(lower chart with reflection at \code{cu}),
and \code{"two"} (two-sided chart), respectively.}
\item{r}{dimension of the resulting linear equation system (highest order of the collocation polynomials).}
\item{qm}{number of quadrature nodes for calculating the collocation definite integrals.}
\item{qm.sigma}{number of quadrature nodes for convoluting the standard deviation uncertainty.}
\item{truncate}{size of truncated tail.}
}
\details{
Essentially, the ARL function \code{sewma.arl} is convoluted with the
distribution of the sample standard deviation.
For details see Jones/Champ/Rigdon (2001) and Knoth (2014?).}
\value{Returns a single value which resembles the ARL.}
\references{
L. A. Jones, C. W. Champ, S. E. Rigdon (2001),
The performance of exponentially weighted moving average charts with estimated parameters,
\emph{Technometrics 43}, 156-167.

S. Knoth (2005),
Accurate ARL computation for EWMA-\eqn{S^2}{S^2} control charts,
\emph{Statistics and Computing 15}, 341-352.

S. Knoth (2006),
Computation of the ARL for CUSUM-\eqn{S^2}{S^2} schemes,
\emph{Computational Statistics & Data Analysis 51}, 499-512.
}
\author{Sven Knoth}
\seealso{
\code{sewma.arl} for zero-state ARL function of EWMA control charts w/o pre run uncertainty.
}
\examples{
## Knoth (2014?)
}
\keyword{ts}
