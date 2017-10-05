\name{scusums.arl}
\alias{scusums.arl}
\title{Compute ARLs of CUSUM-Shewhart control charts (variance charts)}
\description{Computation of the (zero-state) Average Run Length (ARL)
for different types of CUSUM-Shewhart combo control charts (based on the sample variance
\eqn{S^2}) monitoring normal variance.}
\usage{scusums.arl(k, h, cS, sigma, df, hs=0, sided="upper", k2=NULL,
h2=NULL, hs2=0, r=40, qm=30, version=2)}
\arguments{
\item{k}{reference value of the CUSUM control chart.}
\item{h}{decision interval (alarm limit, threshold) of the CUSUM control chart.}
\item{cS}{Shewhart limit.}
\item{sigma}{true standard deviation.}
\item{df}{actual degrees of freedom, corresponds to subgroup size (for known mean it is equal to the subgroup size,
for unknown mean it is equal to subgroup size minus one.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided two-sided CUSUM-\eqn{S^2}{S^2} control charts 
by choosing \code{"upper"} (upper chart), \code{"lower"} (lower chart), and \code{"two"} (two-sided chart), 
respectively. Note that for the two-sided chart the parameters \code{"k2"} and \code{"h2"} have to be set too.}
\item{k2}{In case of a two-sided CUSUM chart for variance the reference value of the lower chart.}
\item{h2}{In case of a two-sided CUSUM chart for variance the decision interval of the lower chart.}
\item{hs2}{In case of a two-sided CUSUM chart for variance the headstart of the lower chart.}
\item{r}{Dimension of the resulting linear equation system (highest order of the collocation
polynomials times number of intervals -- see Knoth 2006).}
\item{qm}{Number of quadrature nodes for calculating the collocation definite integrals.}
\item{version}{Distinguish version numbers (1,2,...). For internal use only.}
}
\details{
\code{scusums.arl} determines the Average Run Length (ARL) by numerically
solving the related ARL integral equation by means of collocation (piecewise Chebyshev polynomials).}
\value{Returns a single value which resembles the ARL.}
\references{
S. Knoth (2006),
Computation of the ARL for CUSUM-\eqn{S^2}{S^2} schemes,
\emph{Computational Statistics & Data Analysis 51}, 499-512.
}
\author{Sven Knoth}
\seealso{
\code{scusum.arl} for zero-state ARL computation of standalone CUSUM control charts for monitoring normal variance.
}
\examples{
## will follow
}
\keyword{ts}
