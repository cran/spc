\name{scusum.crit}
\alias{scusum.crit}
\title{Compute decision intervals of CUSUM control charts (variance charts)}
\description{omputation of the decision intervals (alarm limits)
for different types of CUSUM control charts (based on the sample
variance \eqn{S^2}) monitoring normal variance.}
\usage{scusum.crit(k, L0, sigma, df, hs=0, sided="upper", mode="eq.tails",
k2=NULL, hs2=0, r=40, qm=30)}
\arguments{
\item{k}{reference value of the CUSUM control chart.}
\item{L0}{in-control ARL.}
\item{sigma}{true standard deviation.}
\item{df}{actual degrees of freedom, corresponds to subgroup size (for known mean it is equal
to the subgroup size, for unknown mean it is equal to subgroup size minus one.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided two-sided CUSUM-\eqn{S^2}{S^2} control charts 
by choosing \code{"upper"} (upper chart), \code{"lower"} (lower chart),
and \code{"two"} (two-sided chart), respectively. Note that for the two-sided chart the parameters
\code{"k2"} and \code{"h2"} have to be set too.}
\item{mode}{only deployed for \code{sided}=\code{"two"} -- with \code{"eq.tails"}
two one-sided CUSUM charts (lower and upper) with the same in-control ARL are coupled.
With \code{"unbiased"} a certain unbiasedness of the ARL function is guaranteed (here, both the
lower and the upper control limit are calculated).}
\item{k2}{in case of a two-sided CUSUM chart for variance the reference value of the lower chart.}
\item{hs2}{in case of a two-sided CUSUM chart for variance the headstart of the lower chart.}
\item{r}{Dimension of the resulting linear equation system (highest order of the collocation
polynomials times number of intervals -- see Knoth 2006).}
\item{qm}{Number of quadrature nodes for calculating the collocation definite integrals.}
}
\details{
\code{scusum.crit} ddetermines the decision interval (alarm limit)
for given in-control ARL \code{L0} by applying secant rule and using \code{scusum.arl()}.}
\value{Returns a single value which resembles the decision interval \code{h}.}
\references{
S. Knoth (2005),
Accurate ARL computation for EWMA-\eqn{S^2}{S^2} control charts,
\emph{Statistics and Computing 15}, 341-352.

S. Knoth (2006),
Computation of the ARL for CUSUM-\eqn{S^2}{S^2} schemes,
\emph{Computational Statistics & Data Analysis 51}, 499-512.
}
\author{Sven Knoth}
\seealso{
\code{xcusum.arl} for zero-state ARL computation of CUSUM control charts monitoring normal mean.
}
\examples{
## Knoth (2006)
## compare with Table 1 (p. 507)
k <- 1.46 # sigma1 = 1.5
df <- 1
L0 <- 260.74
h <- scusum.crit(k, L0, 1, df)
h
# original value is 10
}
\keyword{ts}
