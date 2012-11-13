\name{scusum.arl}
\alias{scusum.arl}
\title{Compute ARLs of CUSUM control charts (variance charts)}
\description{Computation of the (zero-state) Average Run Length (ARL)
for different types of CUSUM control charts (based on the sample variance
\eqn{S^2}) monitoring normal variance.}
\usage{scusum.arl(k, h, sigma, df, hs=0, sided="upper", k2=NULL,
h2=NULL, hs2=0, r=40, qm=30, version=2)}
\arguments{
\item{k}{reference value of the CUSUM control chart.}
\item{h}{decision interval (alarm limit, threshold) of the CUSUM control chart.}
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
\code{scusum.arl} determines the Average Run Length (ARL) by numerically
solving the related ARL integral equation by means of collocation (piecewise Chebyshev polynomials).}
\value{Returns a single value which resembles the ARL.}
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
\code{xcusum.arl} for zero-state ARL computation of CUSUM control charts for monitoring normal mean.
}
\examples{
## Knoth (2006)
## compare with Table 1 (p. 507)
k <- 1.46 # sigma1 = 1.5
df <- 1
h <- 10

# original values
# sigma coll63       BE     Hawkins  MC 10^9 (s.e.)
# 1     260.7369  260.7546  261.32  260.7399 (0.0081)
# 1.1    90.1319   90.1389   90.31   90.1319 (0.0027)
# 1.2    43.6867   43.6897   43.75   43.6845 (0.0013)
# 1.3    26.2916   26.2932   26.32   26.2929 (0.0007)
# 1.4    18.1231   18.1239   18.14   18.1235 (0.0005)
# 1.5    13.6268   13.6273   13.64   13.6272 (0.0003)
# 2       5.9904    5.9910    5.99    5.9903 (0.0001)
# replicate the column coll63
sigma <- c(1, 1.1, 1.2, 1.3, 1.4, 1.5, 2)
arl <- rep(NA, length(sigma))
for ( i in 1:length(sigma) )
  arl[i] <- round(scusum.arl(k, h, sigma[i], df, r=63, qm=20, version=2), digits=4)
data.frame(sigma, arl)
}
\keyword{ts}
