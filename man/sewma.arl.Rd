\name{sewma.arl}
\alias{sewma.arl}
\title{Compute ARLs of EWMA control charts (variance charts)}
\description{Computation of the (zero-state) Average Run Length (ARL)
for different types of EWMA control charts (based on the sample variance
\eqn{S^2}) monitoring normal variance.}
\usage{sewma.arl(l,cl,cu,sigma,df,s2.on=TRUE,hs=NULL,sided="upper",r=40,qm=30)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{cl}{lower control limit of the EWMA control chart.}
\item{cu}{upper control limit of the EWMA control chart.}
\item{sigma}{true standard deviation.}
\item{df}{actual degrees of freedom, corresponds to subgroup
size (for known mean it is equal to the subgroup size,
for unknown mean it is equal to subgroup size minus one.}
\item{s2.on}{distinguishes between \eqn{S^2}{S^2} and \eqn{S}{S} chart.}
\item{hs}{so-called headstart (enables fast initial response);
the default (\code{NULL}) yields the expected in-control value of
\eqn{S^2}{S^2} (1) and \eqn{S}{S} (\eqn{c_4}{c_4}), respectively.}
\item{sided}{distinguishes between one- and two-sided 
two-sided EWMA-\eqn{S^2}{S^2} control charts 
by choosing \code{"upper"} (upper chart
without reflection at \code{cl} -- the actual value of
\code{cl} is not used),
\code{"Rupper"} (upper chart with reflection at \code{cl}),
\code{"Rlower"} (lower chart with reflection at \code{cu}),
and \code{"two"} (two-sided chart), 
respectively.} 
\item{r}{dimension of the resulting linear equation system (highest order of
the collocation polynomials).}
\item{qm}{number of quadrature nodes for calculating the collocation
definite integrals.}
}
\details{
\code{sewma.arl} determines the Average Run Length (ARL) by numerically
solving the related ARL integral equation by means of 
collocation (Chebyshev polynomials).}
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
\code{xewma.arl} for zero-state ARL computation of EWMA control charts
for monitoring normal mean.
}
\examples{
## Knoth (2005)
## compare with Table 1 (p. 347): 249.9997
## Monte Carlo with 10^9 replicates: 249.9892 +/- 0.008
l <- .025
df <- 1
cu <- 1 + 1.661865*sqrt(l/(2-l))*sqrt(2/df)
sewma.arl(l,0,cu,1,df)

## ARL values for upper and lower EWMA charts with reflecting barriers
## (reflection at in-control level sigma0 = 1)
## examples from Knoth (2006), Tables 4 and 5

Ssewma.arl <- Vectorize("sewma.arl", "sigma")

## upper chart with reflection at sigma0=1 in Table 4
## original entries are
# sigma   ARL
# 1       100.0
# 1.01    85.3
# 1.02    73.4
# 1.03    63.5
# 1.04    55.4
# 1.05    48.7
# 1.1     27.9
# 1.2     12.9
# 1.3     7.86
# 1.4     5.57
# 1.5     4.30
# 2       2.11

\dontrun{
l <- 0.15
df <- 4
cu <- 1 + 2.4831*sqrt(l/(2-l))*sqrt(2/df)
sigmas <- c(1 + (0:5)/100, 1 + (1:5)/10, 2)
arls <- round(Ssewma.arl(l, 1, cu, sigmas, df, sided="Rupper", r=100), digits=2)
data.frame(sigmas, arls)}

## lower chart with reflection at sigma0=1 in Table 5
## original entries are
# sigma   ARL
# 1       200.04
# 0.9     38.47
# 0.8     14.63
# 0.7     8.65
# 0.6     6.31

\dontrun{
l <- 0.115
df <- 5
cl <- 1 - 2.0613*sqrt(l/(2-l))*sqrt(2/df)
sigmas <- c((10:6)/10)
arls <- round(Ssewma.arl(l, cl, 1, sigmas, df, sided="Rlower", r=100), digits=2)
data.frame(sigmas, arls)}
}
\keyword{ts}
